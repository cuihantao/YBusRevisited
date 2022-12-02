include("makeYbus.jl");


using YBusRevisited: SysConfig, SimConfig


function g_update_turbo!(
    model::LineFull,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    param = model.param
    rhs_g = model.rhs_g
    service = model.service
    aux = model.aux

    (; a1, a2, v1, v2) = algeb
    (; phi, u) = param
    (; a1_rhs, a2_rhs, v1_rhs, v2_rhs) = rhs_g
    (; bh, bhk, gh, ghk, yh, yhk, itap, itap2, cosine, sine) = service
    (;
        itap2_yhyhkconj,
        itap_yhkconj,
        expangle,
        S1_rhs,
        S2_rhs,
        yhconj,
        yhkconj,
        yhyhkconj,
    ) = aux
    a1_ = get_tmp(a1, type_indicator)
    a2_ = get_tmp(a2, type_indicator)
    v1_ = get_tmp(v1, type_indicator)
    v2_ = get_tmp(v2, type_indicator)
    phi_ = get_tmp(phi, type_indicator)
    u_ = get_tmp(u, type_indicator)
    a1_rhs_ = get_tmp(a1_rhs, type_indicator)
    a2_rhs_ = get_tmp(a2_rhs, type_indicator)
    v1_rhs_ = get_tmp(v1_rhs, type_indicator)
    v2_rhs_ = get_tmp(v2_rhs, type_indicator)
    bh_ = get_tmp(bh, type_indicator)
    bhk_ = get_tmp(bhk, type_indicator)
    gh_ = get_tmp(gh, type_indicator)
    ghk_ = get_tmp(ghk, type_indicator)
    itap_ = get_tmp(itap, type_indicator)
    itap2_ = get_tmp(itap2, type_indicator)

    g_update_kernel_turbo!(
        LineFull,
        a1_rhs_,
        a2_rhs_,
        v1_rhs_,
        v2_rhs_,
        a1_,
        a2_,
        bh_,
        bhk_,
        gh_,
        ghk_,
        itap_,
        itap2_,
        phi_,
        u_,
        v1_,
        v2_,
        yhconj,
        yhkconj,
        yhyhkconj,
        expangle,
        S1_rhs,
        S2_rhs,
        itap2_yhyhkconj,
        itap_yhkconj,
    )
    nothing
end


function g_update_kernel_turbo!(
    ::Type{LineFull},
    a1_rhs,
    a2_rhs,
    v1_rhs,
    v2_rhs,
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
    v1,
    v2,
    yhconj,
    yhkconj,
    yhyhkconj,
    expangle,
    S1_rhs,
    S2_rhs,
    itap2_yhyhkconj,
    itap_yhkconj,
)

    @turbo for m in eachindex(u)
        expangle.im[m], expangle.re[m] = sincos(a1[m] - a2[m] - phi[m])
        S1_rhs.re[m] =
            u[m] *
            v1[m] *
            (
                v1[m] * itap2_yhyhkconj.re[m] -
                v2[m] * expangle.re[m] * itap_yhkconj.re[m] +
                v2[m] * expangle.im[m] * itap_yhkconj.im[m]
            )
        S1_rhs.im[m] =
            u[m] *
            v1[m] *
            (
                v1[m] * itap2_yhyhkconj.im[m] -
                v2[m] * expangle.re[m] * itap_yhkconj.im[m] -
                v2[m] * expangle.im[m] * itap_yhkconj.re[m]
            )
        S2_rhs.re[m] =
            u[m] *
            v2[m] *
            (
                v2[m] * yhyhkconj.re[m] -
                v1[m] *
                (itap_yhkconj.re[m] * expangle.re[m] + itap_yhkconj.im[m] * expangle.im[m])
            )
        S2_rhs.im[m] =
            u[m] *
            v2[m] *
            (
                v2[m] * yhyhkconj.im[m] -
                v1[m] * (
                    -itap_yhkconj.re[m] * expangle.im[m] +
                    itap_yhkconj.im[m] * expangle.re[m]
                )
            )
    end

    nothing
end


"""
Collect per-line power injections into per-bus power injections.
"""
function call_collect_inj_S(PQinj, lf::LineFull, addr)
    (; a1_addr, a2_addr, v1_addr, v2_addr) = addr[:Line].algeb
    (; S1_rhs, S2_rhs) = lf.aux

    a1_rhs = S1_rhs.re
    v1_rhs = S1_rhs.im
    a2_rhs = S2_rhs.re
    v2_rhs = S2_rhs.im

    @turbo @. PQinj = 0

    collect_inj_for!(
        PQinj,
        a1_rhs,
        a2_rhs,
        v1_rhs,
        v2_rhs,
        a1_addr,
        a2_addr,
        v1_addr,
        v2_addr,
    )

    nothing
end


@inline function collect_inj_for!(
    PQinj,
    a1_rhs,
    a2_rhs,
    v1_rhs,
    v2_rhs,
    a1_addr,
    a2_addr,
    v1_addr,
    v2_addr,
)

    @inbounds for m in eachindex(a1_rhs)
        PQinj[a1_addr[m]] += a1_rhs[m]
        PQinj[a2_addr[m]] += a2_rhs[m]
        PQinj[v1_addr[m]] += v1_rhs[m]
        PQinj[v2_addr[m]] += v2_rhs[m]
    end
    nothing
end


@views function find_perm_g!(l_nonunique, r_nonunique, rperm, addr)
    rperm[:] .= 0
    l_nonunique .= 0
    r_nonunique .= 0
    idx_nonunique = 1

    for i in eachindex(addr)
        dest = addr[i]
        if rperm[dest] == 0
            rperm[dest] = i
        else
            l_nonunique[idx_nonunique] = dest
            r_nonunique[idx_nonunique] = i
            idx_nonunique += 1
        end
    end
end
