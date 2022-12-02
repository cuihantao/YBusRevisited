"""
Slack Union-type parameters
"""
struct SlackUnionParam
    idx::Vector{String}
    name::Vector{String}
    subidx::Vector{String}
    bus::Vector{String}
    busr::Vector{String}
end

"""
Slack parameter definitions
"""
struct SlackParam{T}
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
    a0::T
    busv0::T
    busa0::T
end

struct SlackState{T} end


SlackState() = SlackState{Nothing}()


struct SlackStateAddr end


struct SlackAlgeb{T}
    q::T
    p::T
    a::T
    v::T
end


struct SlackAlgebAddr{T<:AbstractArray}
    q_addr::T
    p_addr::T
    a_addr::T
    v_addr::T
end


struct SlackRhsF{T} end


SlackRhsF() = SlackRhsF{Nothing}()


struct SlackRhsFAddr end


struct SlackRhsG{T}
    q_rhs::T
    p_rhs::T
    a_rhs::T
    v_rhs::T
end


struct SlackRhsGAddr{T<:AbstractArray}
    q_addr::T
    p_addr::T
    a_addr::T
    v_addr::T
end


struct SlackService{} end


struct SlackDiscrete{T}
    qlim_zi::T
    qlim_zl::T
    qlim_zu::T
    qlim_ql::T
    qlim_qu::T
    plim_zi::T
    plim_zl::T
    plim_zu::T
    plim_ql::T
    plim_qu::T
end


struct SlackConfig
    allow_adjust::Float64
    adjust_lower::Float64
    adjust_upper::Float64
    pv2pq::Float64
    npv2pq::Float64
    min_iter::Float64
    err_tol::Float64
    abs_violation::Float64
    av2pv::Float64
end


struct SlackJacFx{T} end


struct SlackJacFy{T} end


struct SlackJacGx{T} end


struct SlackJacGy{T}
    _gy1::T
    _gy2::T
    _gy3::T
    _gy4::T
    _gy5::T
    _gy6::T
end


"""

    Slack generator.
    
"""
struct Slack{P,S,A,F,G,SV,D}
    n::Int64
    uparam::SlackUnionParam
    param::SlackParam{P}
    state::SlackState{S}
    algeb::SlackAlgeb{A}
    rhs_f::SlackRhsF{F}
    rhs_g::SlackRhsG{G}
    service::SV
    discrete::D
    config::SlackConfig
end


struct SlackAddr
    state::SlackStateAddr
    algeb::SlackAlgebAddr
    rhs_f::SlackRhsFAddr
    rhs_g::SlackRhsGAddr
end


function f_update!(
    model::Slack,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function g_update!(
    model::Slack,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    discrete = model.discrete
    param = model.param
    rhs_g = model.rhs_g

    (; a, p, q, v) = algeb
    (; plim_zi, plim_zl, plim_zu, qlim_zi, qlim_zl, qlim_zu) = discrete
    (; a0, pmax, pmin, qmax, qmin, u, v0) = param
    (; q_rhs, p_rhs, a_rhs, v_rhs) = rhs_g
    a_ = get_tmp(a, type_indicator)
    p_ = get_tmp(p, type_indicator)
    q_ = get_tmp(q, type_indicator)
    v_ = get_tmp(v, type_indicator)
    plim_zi_ = get_tmp(plim_zi, type_indicator)
    plim_zl_ = get_tmp(plim_zl, type_indicator)
    plim_zu_ = get_tmp(plim_zu, type_indicator)
    qlim_zi_ = get_tmp(qlim_zi, type_indicator)
    qlim_zl_ = get_tmp(qlim_zl, type_indicator)
    qlim_zu_ = get_tmp(qlim_zu, type_indicator)
    a0_ = get_tmp(a0, type_indicator)
    pmax_ = get_tmp(pmax, type_indicator)
    pmin_ = get_tmp(pmin, type_indicator)
    qmax_ = get_tmp(qmax, type_indicator)
    qmin_ = get_tmp(qmin, type_indicator)
    u_ = get_tmp(u, type_indicator)
    v0_ = get_tmp(v0, type_indicator)
    q_rhs_ = get_tmp(q_rhs, type_indicator)
    p_rhs_ = get_tmp(p_rhs, type_indicator)
    a_rhs_ = get_tmp(a_rhs, type_indicator)
    v_rhs_ = get_tmp(v_rhs, type_indicator)

    g_update_kernel!(
        Slack,
        q_rhs_,
        p_rhs_,
        a_rhs_,
        v_rhs_,
        a_,
        a0_,
        p_,
        plim_zi_,
        plim_zl_,
        plim_zu_,
        pmax_,
        pmin_,
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
    ::Type{Slack},
    q_rhs,
    p_rhs,
    a_rhs,
    v_rhs,
    a,
    a0,
    p,
    plim_zi,
    plim_zl,
    plim_zu,
    pmax,
    pmin,
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
    @. p_rhs = u * (plim_zi * (-a + a0) + plim_zl * (-p + pmin) + plim_zu * (-p + pmax))
    @. a_rhs = -p * u
    @. v_rhs = -q * u
    nothing
end


@views function collect_f!(dae, model::Slack, addr::SlackAddr, type_ind)
    nothing
end


@views function collect_g!(dae, model::Slack, addr::SlackAddr, type_ind)
    (; q_rhs, p_rhs, a_rhs, v_rhs) = model.rhs_g
    (; q_addr, p_addr, a_addr, v_addr) = addr.rhs_g
    get_tmp(dae.g, type_ind)[q_addr] .+= get_tmp(q_rhs, type_ind)
    get_tmp(dae.g, type_ind)[p_addr] .+= get_tmp(p_rhs, type_ind)
    get_tmp(dae.g, type_ind)[a_addr] .+= get_tmp(a_rhs, type_ind)
    get_tmp(dae.g, type_ind)[v_addr] .+= get_tmp(v_rhs, type_ind)
    nothing
end


@views function deliver_x!(dae, model::Slack, addr::SlackAddr, type_ind)
    nothing
end


@views function deliver_y!(dae, model::Slack, addr::SlackAddr, type_ind)
    (; q, p, a, v) = model.algeb
    (; q_addr, p_addr, a_addr, v_addr) = addr.algeb
    get_tmp(q, type_ind) .= get_tmp(dae.y, type_ind)[q_addr]
    get_tmp(p, type_ind) .= get_tmp(dae.y, type_ind)[p_addr]
    get_tmp(a, type_ind) .= get_tmp(dae.y, type_ind)[a_addr]
    get_tmp(v, type_ind) .= get_tmp(dae.y, type_ind)[v_addr]
    nothing
end


function fx_update!(
    model::Slack,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function fy_update!(
    model::Slack,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gx_update!(
    model::Slack,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gy_update!(
    model::Slack,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    discrete = model.discrete
    jac_gy = model.jac_gy
    param = model.param

    (; plim_zi, plim_zl, plim_zu, qlim_zi, qlim_zl, qlim_zu) = discrete
    (; _gy1, _gy2, _gy3, _gy4, _gy5, _gy6) = jac_gy
    (; u) = param

    gy_update_kernel!(
        Slack,
        _gy1,
        _gy2,
        _gy3,
        _gy4,
        _gy5,
        _gy6,
        plim_zi,
        plim_zl,
        plim_zu,
        qlim_zi,
        qlim_zl,
        qlim_zu,
        u,
    )
    nothing
end


function gy_update_kernel!(
    ::Type{Slack},
    _gy1,
    _gy2,
    _gy3,
    _gy4,
    _gy5,
    _gy6,
    plim_zi,
    plim_zl,
    plim_zu,
    qlim_zi,
    qlim_zl,
    qlim_zu,
    u,
)
    @. _gy1 = u * (-qlim_zl - qlim_zu)
    @. _gy2 = -qlim_zi * u
    @. _gy3 = u * (-plim_zl - plim_zu)
    @. _gy4 = -plim_zi * u
    @. _gy5 = -u
    @. _gy6 = -u
    nothing
end


function l_update_var!(model::Slack; dae_t = 0.0)
    # --- expand structs ---
    algeb = model.algeb
    discrete = model.discrete
    param = model.param

    (; q, p) = algeb
    (;
        qlim_zi,
        qlim_zl,
        qlim_zu,
        qlim_ql,
        qlim_qu,
        plim_zi,
        plim_zl,
        plim_zu,
        plim_ql,
        plim_qu,
    ) = discrete
    (; qmin, qmax, pmin, pmax) = param

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
    l_update_var!(
        Val{:SortedLimiter},
        p,
        pmin,
        pmax,
        plim_zi,
        plim_zl,
        plim_zu,
        plim_ql,
        plim_qu;
        equal = false,
        enable = true,
        no_lower = false,
        no_upper = false,
        sign_lower = 1,
        sign_upper = 1,
    )
    nothing
end


function l_update_eq!(model::Slack; dae_t = 0.0)
    # empty function
    nothing
end
