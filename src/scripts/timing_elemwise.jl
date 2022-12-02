include("base_import.jl");
include("base_loadcase.jl");
include("linefull.jl");
include("makeYbus.jl");

include("g_elemwise.jl");
include("gy_elemwise.jl");


@info "Benchmark: Element-wise method for $case"

bresults_ewise = Vector{BenchmarkTools.Trial}();

# -------------------------------------------------------
# g residual benchmark
# -------------------------------------------------------

@info "Step 1/4: calculate g residual"
g_update_turbo!(lf, sys.config, sys.sim_config, 0);
b = @benchmark g_update_turbo!($lf, $(sys.config), $(sys.sim_config), 0);
push!(bresults_ewise, b);

PQelem = zeros(Float64, 2 * length(V));


@info "Step 2/4: assemble g residual"
call_collect_inj_S(PQelem, lf, (sys.addr));
b = @benchmark call_collect_inj_S($PQelem, $lf, $(sys.addr));
push!(bresults_ewise, b);

# -------------------------------------------------------
# Jacobian gy part element-wise
# -------------------------------------------------------

gy_update_turbo!(lf, sys.config, sys.sim_config, 0)
@info "Step 3/4: calculate gy elements"
b = @benchmark gy_update_turbo!($lf, $sys.config, $sys.sim_config, 0);
push!(bresults_ewise, b);

gy_rows = Vector{Int64}();
gy_cols = similar(gy_rows);
get_row_col_idx!(ss, gy_rows, gy_cols);

gy_vals = lf.jac_gy.mem;

Jelem = sparse(gy_rows, gy_cols, gy_vals, 2 * bus_n, 2 * bus_n);

# --- Separate unique and non-unique indices ---

perm_gy = zeros(Int64, nnz(Jelem));
l_nonunique = zeros(Int64, length(gy_vals) - nnz(Jelem));
r_nonunique = zeros(Int64, length(gy_vals) - nnz(Jelem));
find_perm_gy_vals!(l_nonunique, r_nonunique, perm_gy, Jelem, gy_rows, gy_cols)

assemble_J_elem_2step(Jelem, gy_vals, perm_gy, l_nonunique, r_nonunique)
@info "Step 4/4: assemble gy elements"
b = @benchmark assemble_J_elem_2step($Jelem, $gy_vals, $perm_gy, $l_nonunique, $r_nonunique);
push!(bresults_ewise, b);

@show bresults_ewise
