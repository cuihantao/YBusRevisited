
using SparseArrays

function gy_update_turbo!(
    model::LineFull,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    jac_gy = model.jac_gy
    param = model.param
    service = model.service
    aux = model.aux

    (; a1, a2, v1, v2) = algeb
    (;
        _gy1,
        _gy2,
        _gy3,
        _gy4,
        _gy5,
        _gy6,
        _gy7,
        _gy8,
        _gy9,
        _gy10,
        _gy11,
        _gy12,
        _gy13,
        _gy14,
        _gy15,
        _gy16,
    ) = jac_gy
    (; phi, u) = param
    (;
        bh,
        bhk,
        gh,
        ghk,
        itap,
        itap2,
        v1v2,
        bhksine,
        bhkcosine,
        ghksine,
        ghkcosine,
        itapv1v2,
        itapv1,
        itapv2,
    ) = service
    (; expangle, itap2_yhyhkconj, yhyhkconj) = aux

    (a1, a2, v1, v2) = map(x -> get_tmp(x, type_indicator), (a1, a2, v1, v2))

    # precalculate
    @turbo for m in eachindex(u)
        bhksine[m] = -bhk[m] * expangle.im[m]
        bhkcosine[m] = bhk[m] * expangle.re[m]
        ghksine[m] = -ghk[m] * expangle.im[m]
        ghkcosine[m] = ghk[m] * expangle.re[m]
        itapv1[m] = itap[m] * v1[m]
        itapv2[m] = itap[m] * v2[m]
        itapv1v2[m] = itap[m] * v1[m] * v2[m]
    end

    gy_update_kernel_turbo!(
        LineFull,
        _gy1,
        _gy2,
        _gy3,
        _gy4,
        _gy5,
        _gy6,
        _gy7,
        _gy8,
        _gy9,
        _gy10,
        _gy11,
        _gy12,
        _gy13,
        _gy14,
        _gy15,
        _gy16,
        a1,
        a2,
        bh,
        bhk,
        gh,
        ghk,
        itap,
        itap2,
        phi,
        u,
        itap2_yhyhkconj,
        yhyhkconj,
        v1,
        v2,
        bhksine,
        bhkcosine,
        ghksine,
        ghkcosine,
        itapv1v2,
        itapv1,
        itapv2,
    )
    nothing
end


function gy_update_kernel_turbo!(
    ::Type{LineFull},
    _gy1,
    _gy2,
    _gy3,
    _gy4,
    _gy5,
    _gy6,
    _gy7,
    _gy8,
    _gy9,
    _gy10,
    _gy11,
    _gy12,
    _gy13,
    _gy14,
    _gy15,
    _gy16,
    a1,
    a2,
    bh,
    bhk,
    gh,
    ghk,
    itap,
    itap2,
    phi,
    u,
    itap2_yhyhkconj,
    yhyhkconj,
    v1,
    v2,
    bhksine,
    bhkcosine,
    ghksine,
    ghkcosine,
    itapv1v2,
    itapv1,
    itapv2,
)
    @turbo for m in eachindex(u)
        _gy1[m] = -u[m] * itapv1v2[m] * (bhkcosine[m] + ghksine[m])
        _gy2[m] = -u[m] * itapv1v2[m] * (-bhkcosine[m] - ghksine[m])
        _gy3[m] =
            u[m] *
            (-itapv2[m] * (-bhksine[m] + ghkcosine[m]) + 2 * v1[m] * itap2_yhyhkconj.re[m])
        _gy4[m] = -u[m] * itapv1[m] * (-bhksine[m] + ghkcosine[m])
        _gy5[m] = -u[m] * itapv1v2[m] * (-bhkcosine[m] + ghksine[m])
        _gy6[m] = -u[m] * itapv1v2[m] * (bhkcosine[m] - ghksine[m])
        _gy7[m] = -u[m] * itapv2[m] * (bhksine[m] + ghkcosine[m])
        _gy8[m] =
            u[m] * (-itapv1[m] * (bhksine[m] + ghkcosine[m]) + 2 * v2[m] * yhyhkconj.re[m])
        _gy9[m] = -u[m] * itapv1v2[m] * (-bhksine[m] + ghkcosine[m])
        _gy10[m] = -u[m] * itapv1v2[m] * (bhksine[m] - ghkcosine[m])
        _gy11[m] =
            u[m] *
            (-itapv2[m] * (-bhkcosine[m] - ghksine[m]) + 2 * v1[m] * itap2_yhyhkconj.im[m])
        _gy12[m] = -u[m] * itapv1[m] * (-bhkcosine[m] - ghksine[m])
        _gy13[m] = u[m] * itapv1v2[m] * (bhksine[m] + ghkcosine[m])
        _gy14[m] = u[m] * itapv1v2[m] * (-bhksine[m] - ghkcosine[m])
        _gy15[m] = u[m] * itapv2[m] * (bhkcosine[m] - ghksine[m])
        _gy16[m] =
            u[m] * (itapv1[m] * (bhkcosine[m] - ghksine[m]) + 2 * v2[m] * yhyhkconj.im[m])
    end
    nothing
end


function get_row_col_idx!(ss, rowidx, colidx)

    for i in range(0, 15)
        eqn_name, var_name = ss.Line._jac_eq_var_name("gy", i)

        eqn_addr = getproperty(ss.Line, eqn_name).a .+ 1
        var_addr = getproperty(ss.Line, var_name).a .+ 1
        append!(rowidx, eqn_addr)
        append!(colidx, var_addr)
    end

    nothing

end


function flatten_gy(lf)

    (;
        _gy1,
        _gy2,
        _gy3,
        _gy4,
        _gy5,
        _gy6,
        _gy7,
        _gy8,
        _gy9,
        _gy10,
        _gy11,
        _gy12,
        _gy13,
        _gy14,
        _gy15,
        _gy16,
    ) = lf.jac_gy

    collect(
        Iterators.flatten([
            _gy1,
            _gy2,
            _gy3,
            _gy4,
            _gy5,
            _gy6,
            _gy7,
            _gy8,
            _gy9,
            _gy10,
            _gy11,
            _gy12,
            _gy13,
            _gy14,
            _gy15,
            _gy16,
        ]),
    )

end


@views function update_sparsity_pattern!(pattern, rowidx, colidx, vals)

    pattern.nzval .= 0.0

    @inbounds for idx in eachindex(vals)
        pattern[rowidx[idx], colidx[idx]] += vals[idx]  # values are not unique. need +=
    end
    nothing
end


"""
Find the indices to access `spmat.nzval` in given order.

`perm` stores these indices such that `spmat.nzval[perm]` accesses the values
corresponding to `rowidx` and `colidx`.
"""
function find_nzloc!(perm, spmat, rowidx, colidx)
    for i in eachindex(colidx)
        col = colidx[i]
        col_begin = spmat.colptr[col]
        col_end = spmat.colptr[col+1] - 1
        for k = col_begin:col_end
            if spmat.rowval[k] == rowidx[i]
                perm[i] = k
            end
        end
        if perm[i] == 0
            @error "indices not found"
        end
    end
end


"""
Find permutation matrix to access `gy_vals` in the linear order given by spmat.nzval.

```
spmat.nzval = gy_vals[rperm]

for (l, r) in zip(l_nonunique, r_nonunique)
    spmat.nzval[l] += gy_vals[r]
end
```

"""
function find_perm_gy_vals!(l_nonunique, r_nonunique, rperm, spmat, rowidx, colidx)
    rperm .= 0
    l_nonunique .= 0
    r_nonunique .= 0

    idx_nonunique = 1

    for i in eachindex(colidx)
        col = colidx[i]
        for k in nzrange(spmat, col)
            if spmat.rowval[k] == rowidx[i]
                if rperm[k] == 0
                    rperm[k] = i
                else
                    l_nonunique[idx_nonunique] = k
                    r_nonunique[idx_nonunique] = i
                    idx_nonunique += 1
                end
            end
        end
    end

end


"""
Assemble `gy_vals` into `Jelem` based on `perm`, which is the access order for
`Jelem.nzval`.
"""
function assemble_J_elem!(Jelem, gy_vals, perm)
    Jelem.nzval .= 0
    @turbo for i in eachindex(perm)
        Jelem.nzval[perm[i]] += gy_vals[i]
    end

    # the in-place broadcasting has similar perf. to @turbo
    # @inbounds @views Jelem.nzval[perm] .+= gy_vals
    nothing
end


"""
Assemble `gy_vals` into `Jelem` in two steps.

The first step uses `rperm` to access `gy_vals` and add the resultant array
linearly to `gy_vals.nzval`. @turbo can be used for moderate speed up.

The second step uses `l_nonunique` and `r_nonunique` to add the non-unique
values in place. @turbo must not be to avoid dirty write.
"""
function assemble_J_elem_2step(Jelem, gy_vals, perm_gy, l_nonunique, r_nonunique)
    Jelem.nzval .= 0
    @turbo for i in eachindex(perm_gy)
        Jelem.nzval[i] = gy_vals[perm_gy[i]]
    end

    # must not use @turbo here; results will otherwise be incorrect
    @inbounds for i in eachindex(l_nonunique)
        lidx = l_nonunique[i]
        ridx = r_nonunique[i]
        Jelem.nzval[lidx] += gy_vals[ridx]
    end

    nothing
end
