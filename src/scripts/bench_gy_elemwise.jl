include("base_import.jl");
include("base_loadcase.jl");
include("linefull.jl");
include("gy_elemwise.jl");


#  629.483 μs (4 allocations: 176 bytes)
# gy_update!(lf, sys.config, sys.sim_config, 0)
# @btime gy_update!($lf, $sys.config, $sys.sim_config, 0)


#  431.657 μs (4 allocations: 176 bytes)
# gy_update_optimized1!(lf, sys.config, sys.sim_config, 0)
# @btime gy_update_optimized1!($lf, $sys.config, $sys.sim_config, 0)


# 9241
#  402.111 μs (4 allocations: 176 bytes)
# 70k
#   3.530 ms (0 allocations: 0 bytes)
gy_update_optimized2!(lf, sys.config, sys.sim_config, 0)
@btime gy_update_optimized2!($lf, $sys.config, $sys.sim_config, 0)


# -----------------------------------------
# ---- Assemble into a Jacobian matrix ---
# -----------------------------------------


gy_update_optimized2!(lf, sys.config, sys.sim_config, 0)

rowidx = Vector{Int64}();
colidx = similar(rowidx);

rowidx, colidx = get_row_col_idx!(ss, rowidx, colidx);
gy = flatten_gy(lf);

bus_n = sys.data[:Bus].n;


# -----------------------------------------------------
# Updating gy in place
# very slow here
# -----------------------------------------------------


# 9241
#   8.212 ms (0 allocations: 0 bytes)
pattern = sparse(rowidx, colidx, 0.0, 2 * bus_n, 2 * bus_n);
update_sparsity_pattern!(pattern, rowidx, colidx, gy);
@btime update_sparsity_pattern!($pattern, $rowidx, $colidx, $gy);


# --- This reset time can be neglected ---
#  39.938 μs (0 allocations: 0 bytes)
# @btime $(pattern.nzval) .= 0


# -----------------------------------------------------
# Construct a new Jacobian matrix
# -----------------------------------------------------

# 9241
#   3.754 ms (20 allocations: 6.64 MiB)
# 70k
#   20.461 ms (20 allocations: 39.18 MiB)
@benchmark sparse($gy_rows, $gy_cols, $gy_vals, 2 * bus_n, 2 * bus_n)

Jelem = sparse(rowidx, colidx, gy, 2 * bus_n, 2 * bus_n);


# -----------------------------------------------------
# Sort elements and add to gy inplace
# -----------------------------------------------------

# p = sortperm(colidx);
# rowidx .= rowidx[p];
# colidx .= colidx[p];
# gy .= gy[p];
# @btime begin
# gy .= gy[p];
# end

#   4.562 ms (0 allocations: 0 bytes)
update_sparsity_pattern!(pattern, rowidx, colidx, gy);
@btime update_sparsity_pattern!($pattern, $rowidx, $colidx, $gy);


# --- construct new sparse matrix from sorted indices ---
#  3.868 ms (20 allocations: 6.64 MiB)
@btime sparse($rowidx, $colidx, $gy, 2 * bus_n, 2 * bus_n);

# ---------------------------------------------
# Verify correctness with ANDES
# ---------------------------------------------

ss.PFlow.nr_step()
row_andes = collect(Iterators.flatten(ss.Line.triplets.ijac["gy"]));
col_andes = collect(Iterators.flatten(ss.Line.triplets.jjac["gy"]));
val_andes = collect(Iterators.flatten(ss.Line.triplets.vjac["gy"]));

@assert all(row_andes .== rowidx .- 1) == true
@assert all(col_andes .== colidx .- 1) == true
@assert maximum(abs.(val_andes - gy)) < 1e-10

sparse_gy = sparse(rowidx, colidx, gy, 2 * bus_n, 2 * bus_n);
