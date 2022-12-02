include("base_import.jl");
include("base_loadcase.jl");
include("linefull.jl");
include("makeYbus.jl");
include("gy_csc.jl");


Yx = Ybus.nzval;
Yp = Ybus.colptr;
Yi = Ybus.rowval;

E = U;

bus_n = length(V);
Ibus = StructArray(zeros(ComplexF64, length(V)));
buffer = StructArray(zeros(ComplexF64, length(V)));

dS_dVm = copy(Ybus);
dS_dVa = copy(Ybus);

# @btime dS_dVa = copy($Ybus);


dS_dVm, dS_dVa = dSbus_dV_csc(Ybus, dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yi, V, E);

dS = diagVc * conj(Ybus * diagVn) .+ conj(diagIc) * diagVn;              # dS_dVm
dR = 1im * conj(conj(diagVc) * (diagIc - Ybus * diagVc));                # dS_dVa

droptol!(dS_dVm - dS, 1e-8)
droptol!(dS_dVa - dR, 1e-8)

# -------------------------------------------------------
# Benchmark time to compute complex Jacobians with CSC
#
# 9241
#  584.069 μs (0 allocations: 0 bytes)
#
# 70k
#   3.158 ms (0 allocations: 0 bytes)
#
# --------------------------------------------------------

@benchmark dS_dVm, dS_dVa =
    dSbus_dV_csc($Ybus, $dS_dVm, $dS_dVa, $buffer, $Ibus, $Yx, $Yp, $Yi, $V, $E)
