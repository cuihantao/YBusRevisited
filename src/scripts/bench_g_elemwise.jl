include("base_import.jl");
include("base_loadcase.jl");

include("g_elemwise.jl");


# -------------------------------------
# Benchmark g_update! unoptimized version
#    982.597 μs (0 allocations: 0 bytes)
# -------------------------------------
# g_update!(lf, sys.config, sys.sim_config, 0)
# @btime g_update!($lf, $sys.config, $sys.sim_config, 0)


# -------------------------------------
# Benchmark g_update!  with @turbo only
#   206.242 μs (0 allocations: 0 bytes)
# -------------------------------------
# g_update_turbo!(lf, sys.config, sys.sim_config, 0)
# @btime g_update_turbo!($lf, $sys.config, $sys.sim_config, 0)


# **********************************************************************
# -------------------------------------
# Benchmark g_update! with @turbo and three hand-made optimizations
#  Fastest one so far
#
# 9241
#     176.695 μs (0 allocations: 0 bytes)
# 70k
#   1.335 ms (0 allocations: 0 bytes)
# -------------------------------------
# g_update_turbo_optimized1!(lf, sys.config, sys.sim_config, 0)
# @btime g_update_turbo_optimized1!($lf, $sys.config, $sys.sim_config, 0)


# -------------------------------------
# Benchmark g_update! with optimized version 2;
#
#  optimized2 is slower than `optimized1`
#   more time for memory access
#
# 9241
#   189.046 μs (0 allocations: 0 bytes)
#
# 70k
#    1.586 ms (0 allocations: 0 bytes)
# -------------------------------------

g_update_turbo_optimized2!(lf, sys.config, sys.sim_config, 0)
# @btime g_update_turbo_optimized2!($lf, $sys.config, $sys.sim_config, 0)

a1_rhs_optimized2 = copy(get_tmp(lf.rhs_g.a1_rhs, 0));
a2_rhs_optimized2 = copy(get_tmp(lf.rhs_g.a2_rhs, 0));
v1_rhs_optimized2 = copy(get_tmp(lf.rhs_g.v1_rhs, 0));
v2_rhs_optimized2 = copy(get_tmp(lf.rhs_g.v2_rhs, 0));


# -------------------------------------
# Benchmark g_update!  with complex numbers
#
# 9241
#    197.741 μs (0 allocations: 0 bytes) without StructArray storage
#    139.845 μs (0 allocations: 0 bytes)
#    136.109 μs (0 allocations: 0 bytes)  without `u *`
#
# 70k
# -------------------------------------
g_update_turbo_optimized3!(lf, sys.config, sys.sim_config, 0)

S1_rhs = lf.aux.S1_rhs;
S2_rhs = lf.aux.S2_rhs;

a1_rhs_complex = S1_rhs.re;
v1_rhs_complex = S1_rhs.im;
a2_rhs_complex = S2_rhs.re;
v2_rhs_complex = S2_rhs.im;

maximum(abs.(a1_rhs_optimized2 - a1_rhs_complex)) < 1e-10
maximum(abs.(v1_rhs_optimized2 - v1_rhs_complex)) < 1e-10
maximum(abs.(a2_rhs_optimized2 - a2_rhs_complex)) < 1e-10
maximum(abs.(v2_rhs_optimized2 - v2_rhs_complex)) < 1e-10

@btime g_update_turbo_optimized3!($lf, $sys.config, $sys.sim_config, 0);


# --------------------------------------------------
# sum up injections from all lines at buses
#
# Summation time is NOT negligible
#
# 9241
#   77.318 μs (5 allocations: 224 bytes)
#
#  70k
#   465.132 μs (5 allocations: 224 bytes)
# --------------------------------------------------


PQinj = zeros(Float64, 2 * length(V));
call_collect_inj(PQinj, lf, sys.addr);
@btime call_collect_inj($PQinj, $lf, $sys.addr);


include("bench_g_Ybus.jl")

maximum(abs.(PQYbus - PQinj))

# ---------------------------------------------------
# Test if an incidence matrix will make it faster
#
# No. Multiplying by an incidence matrix is a lot slower.
# A single matrix-vector product can take
#    132.294 μs (2 allocations: 546.92 KiB)
#

line_n = lf.n
incidence = spdiagm(bus_n, line_n, ones(bus_n));

a1_rhs = get_tmp(lf.rhs_g.a1_rhs, 0);
a1_addr = get_tmp(sys.addr[:Line].rhs_g.a1_addr, 0);

incidence * a1_rhs

@btime begin
    $incidence * $a1_rhs
end
