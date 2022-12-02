"""
PV Union-type parameters
"""
struct PVUnionParam
    idx::Vector{String}
    name::Vector{String}
    subidx::Vector{String}
    bus::Vector{String}
    busr::Vector{String}
end

"""
PV parameter definitions
"""
struct PVParam{T}
    u::T
    Sn::T
    Vn::T
    p0::T
    q0::T
    pmax::T
    pmin::T
    qmax::T
    qmin::T
    v0::T
    vmax::T
    vmin::T
    ra::T
    xs::T
    busv0::T
end

struct PVState{T} end


PVState() = PVState{Nothing}()


struct PVStateAddr end


struct PVAlgeb{T}
    q::T
    a::T
    v::T
end


struct PVAlgebAddr{T<:AbstractArray}
    q_addr::T
    a_addr::T
    v_addr::T
end


struct PVRhsF{T} end


PVRhsF() = PVRhsF{Nothing}()


struct PVRhsFAddr end


struct PVRhsG{T}
    q_rhs::T
    a_rhs::T
    v_rhs::T
end


struct PVRhsGAddr{T<:AbstractArray}
    q_addr::T
    a_addr::T
    v_addr::T
end


struct PVService{T}
    p::T
end


struct PVDiscrete{T}
    qlim_zi::T
    qlim_zl::T
    qlim_zu::T
    qlim_ql::T
    qlim_qu::T
end


struct PVConfig
    allow_adjust::Float64
    adjust_lower::Float64
    adjust_upper::Float64
    pv2pq::Float64
    npv2pq::Float64
    min_iter::Float64
    err_tol::Float64
    abs_violation::Float64
end


struct PVJacFx{T} end


struct PVJacFy{T} end


struct PVJacGx{T} end


struct PVJacGy{T}
    _gy1::T
    _gy2::T
    _gy3::T
end


"""

    Static PV generator with reactive power limit checking
    and PV-to-PQ conversion.

    `pv2pq = 1` turns on the conversion.
    It starts  from iteration `min_iter` or when the convergence
    error drops below `err_tol`.

    The PV-to-PQ conversion first ranks the reactive violations.
    A maximum number of `npv2pq` PVs above the upper limit, and
    a maximum of `npv2pq` PVs below the lower limit will be
    converted to PQ, which sets the reactive power to `pmax` or
    `pmin`.

    If `pv2pq` is `1` (enabled) and `npv2pq` is `0`, heuristics
    will be used to determine the number of PVs to be converted
    for each iteration.
    
"""
struct PV{P,S,A,F,G,SV,D}
    n::Int64
    uparam::PVUnionParam
    param::PVParam{P}
    state::PVState{S}
    algeb::PVAlgeb{A}
    rhs_f::PVRhsF{F}
    rhs_g::PVRhsG{G}
    service::SV
    discrete::D
    config::PVConfig
end


struct PVAddr
    state::PVStateAddr
    algeb::PVAlgebAddr
    rhs_f::PVRhsFAddr
    rhs_g::PVRhsGAddr
end


function f_update!(model::PV, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function g_update!(model::PV, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    algeb = model.algeb
    discrete = model.discrete
    param = model.param
    rhs_g = model.rhs_g
    service = model.service

    (; q, v) = algeb
    (; qlim_zi, qlim_zl, qlim_zu) = discrete
    (; qmax, qmin, u, v0) = param
    (; q_rhs, a_rhs, v_rhs) = rhs_g
    (; p) = service
    q_ = get_tmp(q, type_indicator)
    v_ = get_tmp(v, type_indicator)
    qlim_zi_ = get_tmp(qlim_zi, type_indicator)
    qlim_zl_ = get_tmp(qlim_zl, type_indicator)
    qlim_zu_ = get_tmp(qlim_zu, type_indicator)
    qmax_ = get_tmp(qmax, type_indicator)
    qmin_ = get_tmp(qmin, type_indicator)
    u_ = get_tmp(u, type_indicator)
    v0_ = get_tmp(v0, type_indicator)
    q_rhs_ = get_tmp(q_rhs, type_indicator)
    a_rhs_ = get_tmp(a_rhs, type_indicator)
    v_rhs_ = get_tmp(v_rhs, type_indicator)
    p_ = get_tmp(p, type_indicator)

    g_update_kernel!(
        PV,
        q_rhs_,
        a_rhs_,
        v_rhs_,
        p_,
        q_,
        qlim_zi_,
        qlim_zl_,
        qlim_zu_,
        qmax_,
        qmin_,
        u_,
        v_,
        v0_,
    )
    nothing
end


function g_update_kernel!(
    ::Type{PV},
    q_rhs,
    a_rhs,
    v_rhs,
    p,
    q,
    qlim_zi,
    qlim_zl,
    qlim_zu,
    qmax,
    qmin,
    u,
    v,
    v0,
)
    @. q_rhs = u * (qlim_zi * (-v + v0) + qlim_zl * (-q + qmin) + qlim_zu * (-q + qmax))
    @. a_rhs = -p * u
    @. v_rhs = -q * u
    nothing
end


@views function collect_f!(dae, model::PV, addr::PVAddr, type_ind)
    nothing
end


@views function collect_g!(dae, model::PV, addr::PVAddr, type_ind)
    (; q_rhs, a_rhs, v_rhs) = model.rhs_g
    (; q_addr, a_addr, v_addr) = addr.rhs_g
    get_tmp(dae.g, type_ind)[q_addr] .+= get_tmp(q_rhs, type_ind)
    get_tmp(dae.g, type_ind)[a_addr] .+= get_tmp(a_rhs, type_ind)
    get_tmp(dae.g, type_ind)[v_addr] .+= get_tmp(v_rhs, type_ind)
    nothing
end


@views function deliver_x!(dae, model::PV, addr::PVAddr, type_ind)
    nothing
end


@views function deliver_y!(dae, model::PV, addr::PVAddr, type_ind)
    (; q, a, v) = model.algeb
    (; q_addr, a_addr, v_addr) = addr.algeb
    get_tmp(q, type_ind) .= get_tmp(dae.y, type_ind)[q_addr]
    get_tmp(a, type_ind) .= get_tmp(dae.y, type_ind)[a_addr]
    get_tmp(v, type_ind) .= get_tmp(dae.y, type_ind)[v_addr]
    nothing
end


function fx_update!(model::PV, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function fy_update!(model::PV, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function gx_update!(model::PV, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function gy_update!(model::PV, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    discrete = model.discrete
    jac_gy = model.jac_gy
    param = model.param

    (; qlim_zi, qlim_zl, qlim_zu) = discrete
    (; _gy1, _gy2, _gy3) = jac_gy
    (; u) = param

    gy_update_kernel!(PV, _gy1, _gy2, _gy3, qlim_zi, qlim_zl, qlim_zu, u)
    nothing
end


function gy_update_kernel!(::Type{PV}, _gy1, _gy2, _gy3, qlim_zi, qlim_zl, qlim_zu, u)
    @. _gy1 = u * (-qlim_zl - qlim_zu)
    @. _gy2 = -qlim_zi * u
    @. _gy3 = -u
    nothing
end


function l_update_var!(model::PV; dae_t = 0.0)
    # --- expand structs ---
    algeb = model.algeb
    discrete = model.discrete
    param = model.param

    (; q) = algeb
    (; qlim_zi, qlim_zl, qlim_zu, qlim_ql, qlim_qu) = discrete
    (; qmin, qmax) = param

    # --- begin `l_update_var` calls ---
    l_update_var!(
        Val{:SortedLimiter},
        q,
        qmin,
        qmax,
        qlim_zi,
        qlim_zl,
        qlim_zu,
        qlim_ql,
        qlim_qu;
        equal = false,
        enable = true,
        no_lower = false,
        no_upper = false,
        sign_lower = 1,
        sign_upper = 1,
    )
    nothing
end


function l_update_eq!(model::PV; dae_t = 0.0)
    # empty function
    nothing
end
