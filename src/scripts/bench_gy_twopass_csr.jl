
include("base_import.jl");
include("base_loadcase.jl");
include("linefull.jl");
include("makeYbus.jl");

using SparseMatricesCSR

# Reimplementation of
# https://github.com/SanPen/GridCal/blob/aaa4ad6cec0c637f67194e98cf0b615f16b6711b/src/GridCal/Engine/Simulations/PowerFlow/NumericalMethods/derivatives.py#L114-L123


"""
partial derivatives of power injection w.r.t. voltage.
:param Yx: Ybus data in CSC format
:param Yp: Ybus colptr in CSC format
:param Yj: Ybus rowval in CSC format
:param V: Voltage vector
:param E: Normalized voltage vector
:return: dS_dVm, dS_dVa data in CSR format, index pointer and indices are the same as the ones from Ybus
"""
function dSbus_dV_csr(dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yj, V, E)

    # `r`: row index
    for r in range(1, length(Yp) - 1)
        for k in range(Yp[r], Yp[r+1] - 1)
            # `Yj[k]`: column index

            # Ibus = Ybus * V
            buffer[r] += Yx[k] * V[Yj[k]]

            # Ybus * diag(Vnorm)
            dS_dVm.nzval[k] *= E[Yj[k]]

            # Ybus * diag(V)
            dS_dVa.nzval[k] *= V[Yj[k]]
        end

        Ibus[r] += buffer[r]

        # conj(diagIbus) * diagVnorm
        buffer[r] = conj(buffer[r]) * E[r]
    end

    for r in range(1, length(Yp) - 1)
        for k in range(Yp[r], Yp[r+1] - 1)
            # diag(V) * conj(Ybus * diagVnorm)
            dS_dVm.nzval[k] = conj(dS_dVm.nzval[k]) * V[r]

            if r == Yj[k]
                # diagonal elements
                dS_dVa.nzval[k] -= Ibus[r]
                dS_dVm.nzval[k] += buffer[r]
            end

            # 1j * diagV * conj(diagIbus - Ybus * diagV)
            dS_dVa.nzval[k] = conj(-dS_dVa.nzval[k]) * (1im * V[r])
        end
    end

    return dS_dVm, dS_dVa

end


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


# ---
# NOTE: need to transpose twice to get the same matrix in CSR
# ---
YbusT = sparse(transpose(Ybus));
Ybus_csr = SparseMatrixCSR(transpose(YbusT));

Yx = Ybus_csr.nzval;
Yp = Ybus_csr.rowptr;
Yj = Ybus_csr.colval;

E = U;
bus_n = length(V);
Ibus = zeros(ComplexF64, length(V));
buffer = zeros(ComplexF64, length(V));

dS_dVm = copy(Ybus_csr);
dS_dVa = copy(Ybus_csr);


dS_dVm, dS_dVa = dSbus_dV_csr(dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yj, V, E);

dS = diagVc * conj(Ybus * diagVn) .+ conj(diagIc) * diagVn;              # dS_dVm
dR = 1im * conj(conj(diagVc) * (diagIc - Ybus * diagVc));                # dS_dVa

dS_diff = dS_dVm - dS
dR_diff = dS_dVa - dR


# 9241
#   622.920 μs (0 allocations: 0 bytes)
# 70k
#     3.052 ms (0 allocations: 0 bytes)
@btime dSbus_dV_csr($dS_dVm, $dS_dVa, $buffer, $Ibus, $Yx, $Yp, $Yj, $V, $E);


# --------------------------------------------
# Unpack dS_dVm and dS_dVa in 2n*2n sparse matrix with Floats
#
# 9214
#    3.758 ms (283 allocations: 24.38 MiB)
# 70k
#   25.357 ms (285 allocations: 135.92 MiB)
#
# --------------------------------------------


# concat by rows first then cols
# This is slower for CSR matrices
# 9241
#    3.842 ms (283 allocations: 24.38 MiB)

# function to_full_J(dS, dR)
#     sparse_vcat(sparse_hcat(real(dR), real(dS)),
#                 sparse_hcat(imag(dR), imag(dS))
#     )
# end
