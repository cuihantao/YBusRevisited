
# Reimplementation of
# https://github.com/SanPen/GridCal/blob/aaa4ad6cec0c637f67194e98cf0b615f16b6711b/src/GridCal/Engine/Simulations/PowerFlow/NumericalMethods/derivatives.py#L114-L123


@inline function set_zero!(num::Vector{Complex{T}}) where {T}
    @inbounds for m in eachindex(num)
        num[m] = 0
    end
    nothing
end

@inline function set_zero!(num::StructArray{T}) where {T}
    @turbo num.re .= 0
    @turbo num.im .= 0
    nothing
end


"""
Compute the power injection derivatives w.r.t the voltage module and angle
:param Yx: data of Ybus in CSC format
:param Yp: indptr of Ybus in CSC format
:param Yi: indices of Ybus in CSC format
:param V: Voltages vector
:return: dS_dVm, dS_dVa data ordered in the CSC format to match the indices of Ybus
"""

"""
The matrix operations that this is performing are:
diagV = diags(V)
diagE = diags(V / np.abs(V))
Ibus = Ybus * V
diagIbus = diags(Ibus)
dSbus_dVa = 1j * diagV * np.conj(diagIbus - Ybus * diagV)
dSbus_dVm = diagV * np.conj(Ybus * diagE) + np.conj(diagIbus) * diagE
"""
function dSbus_dV_csc(Ybus, dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yi, V, E)

    set_zero!(Ibus)
    set_zero!(buffer)

    dS_dVm.nzval .= Ybus.nzval
    dS_dVa.nzval .= Ybus.nzval

    for j in range(1, length(Yp) - 1)  # for each column ...
        for k in range(Yp[j], Yp[j+1] - 1)  # for each row ...
            # Ibus = Ybus * V
            Ibus[Yi[k]] += Yx[k] * V[j]  # Yx[k] -> Y(i,j)

            # Ybus * diagE
            dS_dVm.nzval[k] *= E[j]

            # Ybus * diag(V)
            dS_dVa.nzval[k] *= V[j]
        end
    end

    for j in range(1, length(Yp) - 1)  # for each column ...

        buffer[j] = conj(Ibus[j]) * E[j]
        for k in range(Yp[j], Yp[j+1] - 1)  # for each row ...

            # row index
            i = Yi[k]

            # diag(V) * conj(Ybus * diagE)
            dS_dVm.nzval[k] = V[i] * conj(dS_dVm.nzval[k])

            if j == i
                # diagonal elements
                dS_dVa.nzval[k] -= Ibus[j]
                dS_dVm.nzval[k] += buffer[j]
            end

            # 1j * diagV * conj(diagIbus - Ybus * diagV)
            dS_dVa.nzval[k] = conj(-dS_dVa.nzval[k]) * (1im * V[i])
        end
    end


    return dS_dVm, dS_dVa
end


# concatenate by cols first then rows
# 9241
#    3.458 ms (283 allocations: 23.68 MiB)
function dSdR_to_J(dS, dR)
    # NOTE: this method drops all zeros
    sparse_hcat(sparse_vcat(real(dR), imag(dR)), sparse_vcat(real(dS), imag(dS)))

end

"""
Assemble dS and dR into J while preserving all zeros in the sparsity pattern.
"""
function dSdR_to_J_keep_zeros(dS, dR)
    nrows, ncols = size(dR)
    rowidx, colidx, values_dR = findnz(dR)
    _, _, values_dS = findnz(dS)

    # column-major construction
    rowidx = [rowidx; rowidx .+ nrows; rowidx; rowidx .+ nrows]
    colidx = [colidx; colidx; colidx .+ ncols; colidx .+ ncols]
    values = [real(values_dR); imag(values_dR); real(values_dS); imag(values_dS)]

    sparse(rowidx, colidx, values, 2 * nrows, 2 * ncols)

end


# dSdR = to_full_J(dS, dR);
# @btime to_full_J($dS, $dR);

function inc_perm_idx_col!(perm, colidx, Ybus)
    for m in eachindex(colidx)
        col = colidx[m]
        nrows = Ybus.colptr[col+1] - Ybus.colptr[col]
        perm[m] += nrows
    end
end

function assemble_J_Ybus(Jshape, dS_dVa, dS_dVm, perm_J1, perm_J2, perm_J3, perm_J4)
    Jshape.nzval .= 0

    # Jshape.nzval[perm_J1] .= real(dS_dVa.nzval);
    # Jshape.nzval[perm_J2] .= imag(dS_dVa.nzval);
    # Jshape.nzval[perm_J3] .= real(dS_dVm.nzval);
    # Jshape.nzval[perm_J4] .= imag(dS_dVm.nzval);

    # for (i, m) in enumerate(perm_J1)
    #     Jshape.nzval[m] = real(dS_dVa.nzval[i]);
    # end
    # for (i, m) in enumerate(perm_J2)
    #     Jshape.nzval[m] = imag(dS_dVa.nzval[i]);
    # end
    # for (i, m) in enumerate(perm_J3)
    #     Jshape.nzval[m] = real(dS_dVm.nzval[i]);
    # end
    # for (i, m) in enumerate(perm_J4)
    #     Jshape.nzval[m] = imag(dS_dVm.nzval[i]);
    # end

    # similar perf. with or without @inbounds

    for i in eachindex(perm_J1)
        Jshape.nzval[perm_J1[i]] = real(dS_dVa.nzval[i])
        Jshape.nzval[perm_J2[i]] = imag(dS_dVa.nzval[i])
        Jshape.nzval[perm_J3[i]] = real(dS_dVm.nzval[i])
        Jshape.nzval[perm_J4[i]] = imag(dS_dVm.nzval[i])
    end

    nothing
end
