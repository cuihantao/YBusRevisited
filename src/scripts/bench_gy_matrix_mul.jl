
include("base_import.jl");
include("base_loadcase.jl");
include("linefull.jl");
include("makeYbus.jl");

#=
n1 = Bus.a;
U = exp(i*DAE.y(n1));
V = DAE.y(Bus.v).*U;
I = a.Y*V;

diagVc = sparse(n1,n1,V,nb,nb);
diagVn = sparse(n1,n1,U,nb,nb);
diagIc = sparse(n1,n1,I,nb,nb);
dS = diagVc * conj(a.Y * diagVn) + conj(diagIc) * diagVn;
dR = conj(diagVc) * (diagIc - a.Y * diagVc);

[h,k,s] = find([imag(dR),real(dS);real(dR),imag(dS)]);
Gy = sparse(h,k,s,DAE.m,DAE.m);
=#

bus = sys.data[:Bus];

bus_n = bus.n;
bus_v = get_tmp(bus.algeb.v, 0);
bus_a = get_tmp(bus.algeb.a, 0);

U = exp.(1im * bus_a);
V = U .* bus_v;
Ic = Ybus * V;

diagVc = spdiagm(bus_n, bus_n, V);
diagVn = spdiagm(bus_n, bus_n, U);
diagIc = spdiagm(bus_n, bus_n, Ic);

# (diagVc, diagVn, diagIc) = map(x->convert(Diagonal, x), (diagVc, diagVn, diagIc));

# --------------------------------------------------------
# Benchmark `dS`
# --------------------------------------------------------

dS = diagVc * conj(Ybus * diagVn) .+ conj(diagIc) * diagVn;  # dS_dVm
dR = conj(diagVc) * (diagIc - Ybus * diagVc);                # dS_dVa


function calc_dS!(dS, diagVc, Ybus, diagVn, diagIc)
    dS .= diagVc * conj(Ybus * diagVn) .+ conj(diagIc) * diagVn
end


calc_dS!(dS, diagVc, Ybus, diagVn, diagIc);
# 9241
#   1.473 ms (33 allocations: 4.22 MiB)
# 70k
#   9.240 ms (36 allocations: 28.37 MiB)
@btime calc_dS!($dS, $diagVc, $Ybus, $diagVn, $diagIc);


out1 = spzeros(Complex, bus_n, bus_n);
out2 = spzeros(Complex, bus_n, bus_n);
out3 = spzeros(Complex, bus_n, bus_n);

"""
This function preallocates memory for intermediate variables but is considerably
slower.

Even without the output line,
  4.438 ms (0 allocations: 0 bytes)

"""
function calc_dS_noalloc!(out1, out2, out3, dS, diagVc, Ybus, diagVn, diagIc)

    # ybus = copy(Ybus)
    ybus = Ybus
    rmul!(ybus, diagVn)
    conj!(ybus)

    lmul!(diagVc, ybus)

    conj!(diagIc)
    lmul!(diagIc, diagVn)

    # dS .= out2 + out3
    nothing

end

calc_dS_noalloc!(out1, out2, out3, dS, diagVc, Ybus, diagVn, diagIc)
@btime calc_dS_noalloc!($out1, $out2, $out3, $dS, $diagVc, $Ybus, $diagVn, $diagIc)
@code_typed calc_dS_noalloc!(out1, out2, out3, dS, diagVc, Ybus, diagVn, diagIc)


# --------------------------------------------------------
# Benchmark `dR`
# --------------------------------------------------------

#  1.215 ms (30 allocations: 3.90 MiB)
@btime dR = conj(diagVc) * (diagIc - Ybus * diagVc);


# --------------------------------------------------------
# CUDA version to calculate `dS`
# --------------------------------------------------------

using CUDA
using CUDA.CUSPARSE


Ybus_d = cu(Ybus);
U_d = cu(U);
V_d = cu(V);
I_d = cu(I);

diagVc_d = cu(sparse(diagVc));
diagVn_d = cu(sparse(diagVn));
diagIc_d = cu(sparse(diagIc));


CUDA.allowscalar(false)

(Ybus_d, diagVc_d, diagVn_d, diagIc_d) =
    map(x -> convert(CuSparseMatrixCSR, x), (Ybus_d, diagVc_d, diagVn_d, diagIc_d));


out1 = CUSPARSE.CuSparseMatrixCSR(spzeros(eltype(Ybus_d), size(Ybus_d, 1), size(Ybus_d, 2)));
out2 = CUSPARSE.CuSparseMatrixCSR(spzeros(eltype(Ybus_d), size(Ybus_d, 1), size(Ybus_d, 2)));
out3 = CUSPARSE.CuSparseMatrixCSR(spzeros(eltype(Ybus_d), size(Ybus_d, 1), size(Ybus_d, 2)));
dS_d = CUSPARSE.CuSparseMatrixCSR(spzeros(eltype(Ybus_d), size(Ybus_d, 1), size(Ybus_d, 2)));


function csrmm!(
    C::CUSPARSE.CuSparseMatrixCSR,
    A::CUSPARSE.CuSparseMatrixCSR,
    B::CUSPARSE.CuSparseMatrixCSR,
)
    CUSPARSE.mm!('N', 'N', one(eltype(A)), A, B, zero(eltype(A)), C, 'O')
    nothing
end


function calc_dS_cuda!(out1, out2, out3, dS_d, diagVc_d, Ybus_d, diagVn_d, diagIc_d)
    csrmm!(out1, Ybus_d, diagVn_d)
    csrmm!(out2, diagVc_d, out1)
    csrmm!(out3, conj(diagIc_d), diagVn_d)
    dS_d = out2 + out3
    nothing
end

calc_dS_cuda!(out1, out2, out3, dS_d, diagVc_d, Ybus_d, diagVn_d, diagIc_d)

#  598.945 μs (435 allocations: 13.56 KiB)
@btime calc_dS_cuda!(out1, out2, out3, dS_d, diagVc_d, Ybus_d, diagVn_d, diagIc_d)
