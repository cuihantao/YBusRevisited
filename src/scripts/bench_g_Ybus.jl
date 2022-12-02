include("g_Ybus.jl")

PQYbus = zeros(Float64, 2 * length(V));
S = similar(V);
bus_n = length(V);


# -----------------------------------------------------
# Ybus method to calculate power injections
# -----------------------------------------------------

g_update_Ybus!(PQYbus, S, V, Ybus, bus_n);

@btime g_update_Ybus!($PQYbus, $S, $V, $Ybus, $bus_n)


using SuiteSparseGraphBLAS


#------------------------------------------------
# GraphBLAS
#
#  case9241
#   203.212 μs (35 allocations: 1.31 KiB)
#
# julia -t 2
#   199.618 μs (35 allocations: 1.31 KiB)  -- # of threads is set to maximum
#------------------------------------------------

Ybus_GB = GBMatrix(Ybus);
V_GB = GBMatrix(V);
S_GB = GBMatrix{Float64}(bus_n, 1, fill = 1.0);
PQYbus_GB = Vector{Float64}(undef, 2 * bus_n);

ret = g_update_Ybus_GB!(PQYbus_GB, S_GB, V_GB, Ybus_GB, bus_n);
@btime g_update_Ybus_GB!($PQYbus_GB, $S_GB, $V_GB, $Ybus_GB, $bus_n);
