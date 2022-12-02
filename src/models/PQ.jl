"""
PQ Union-type parameters
"""
struct PQUnionParam
    idx::Vector{String}
    name::Vector{String}
    bus::Vector{String}
    owner::Vector{String}
end

"""
PQ parameter definitions
"""
struct PQParam{T}
    u::T
    Vn::T
    p0::T
    q0::T
    vmax::T
    vmin::T
end

struct PQState{T} end


PQState() = PQState{Nothing}()


struct PQStateAddr end


struct PQAlgeb{T}
    a::T
    v::T
end


struct PQAlgebAddr{T<:AbstractArray}
    a_addr::T
    v_addr::T
end


struct PQRhsF{T} end


PQRhsF() = PQRhsF{Nothing}()


struct PQRhsFAddr end


struct PQRhsG{T}
    a_rhs::T
    v_rhs::T
end


struct PQRhsGAddr{T<:AbstractArray}
    a_addr::T
    v_addr::T
end


struct PQService{T}
    Rub::T
    Xub::T
    Rlb::T
    Xlb::T
    Ppf::T
    Qpf::T
    Req::T
    Xeq::T
    Ipeq::T
    Iqeq::T
    v0::T
    a0::T
end


struct PQDiscrete{T}
    vcmp_zi::T
    vcmp_zl::T
    vcmp_zu::T
end


struct PQConfig
    allow_adjust::Float64
    adjust_lower::Float64
    adjust_upper::Float64
    pq2z::Float64
    p2p::Float64
    p2i::Float64
    p2z::Float64
    q2q::Float64
    q2i::Float64
    q2z::Float64
end


struct PQJacFx{T} end


struct PQJacFy{T} end


struct PQJacGx{T} end


struct PQJacGy{T}
    _gy1::T
    _gy2::T
end


"""

    PQ load model.

    Implements an automatic pq2z conversion during power flow when the voltage
    is outside [vmin, vmax]. The conversion can be turned off by setting `pq2z`
    to 0 in the Config file.

    Before time-domain simulation, PQ load will be converted to impedance,
    current source, and power source based on the weights in the Config file.

    Weights (p2p, p2i, p2z) corresponds to the weights for constant power,
    constant current and constant impedance. p2p, p2i and p2z must be in decimal
    numbers and sum up exactly to 1. The same rule applies to (q2q, q2i, q2z).

    To alter the PQ load in terms of power during simulation, one needs to set
    the conversion weights to preserve the constant power portion. For example,
    the PQ can remain as constant power load by setting

    .. code-block :: python

        ss.PQ.config.p2p = 1.0 ss.PQ.config.p2i = 0 ss.PQ.config.p2z = 0

        ss.PQ.config.q2q = 1.0 ss.PQ.config.q2i = 0 ss.PQ.config.q2z = 0

    Then, the constant power portion can be altered by changing the ``Ppf`` and
    ``Qpf`` constants for active power and reactive power.

    The equivalent constant current components are in constants ``Ipeq`` and
    ``Iqeq`` for active and reactive current, and the equivalent impedances are
    in ``Req`` and ``Xeq``.
    
"""
struct PQ{P,S,A,F,G,SV,D}
    n::Int64
    uparam::PQUnionParam
    param::PQParam{P}
    state::PQState{S}
    algeb::PQAlgeb{A}
    rhs_f::PQRhsF{F}
    rhs_g::PQRhsG{G}
    service::SV
    discrete::D
    config::PQConfig
end


struct PQAddr
    state::PQStateAddr
    algeb::PQAlgebAddr
    rhs_f::PQRhsFAddr
    rhs_g::PQRhsGAddr
end


function f_update!(model::PQ, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function g_update!(model::PQ, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    algeb = model.algeb
    config = model.config
    discrete = model.discrete
    param = model.param
    rhs_g = model.rhs_g
    service = model.service

    (; dae_t) = sim_config
    (; v) = algeb
    (; p2i, p2p, p2z, q2i, q2q, q2z) = config
    (; vcmp_zi, vcmp_zl, vcmp_zu) = discrete
    (; p0, q0, u) = param
    (; a_rhs, v_rhs) = rhs_g
    (; Ipeq, Iqeq, Ppf, Qpf, Req, Rlb, Rub, Xeq, Xlb, Xub) = service
    v_ = get_tmp(v, type_indicator)
    p2i_ = get_tmp(p2i, type_indicator)
    p2p_ = get_tmp(p2p, type_indicator)
    p2z_ = get_tmp(p2z, type_indicator)
    q2i_ = get_tmp(q2i, type_indicator)
    q2q_ = get_tmp(q2q, type_indicator)
    q2z_ = get_tmp(q2z, type_indicator)
    vcmp_zi_ = get_tmp(vcmp_zi, type_indicator)
    vcmp_zl_ = get_tmp(vcmp_zl, type_indicator)
    vcmp_zu_ = get_tmp(vcmp_zu, type_indicator)
    p0_ = get_tmp(p0, type_indicator)
    q0_ = get_tmp(q0, type_indicator)
    u_ = get_tmp(u, type_indicator)
    a_rhs_ = get_tmp(a_rhs, type_indicator)
    v_rhs_ = get_tmp(v_rhs, type_indicator)
    Ipeq_ = get_tmp(Ipeq, type_indicator)
    Iqeq_ = get_tmp(Iqeq, type_indicator)
    Ppf_ = get_tmp(Ppf, type_indicator)
    Qpf_ = get_tmp(Qpf, type_indicator)
    Req_ = get_tmp(Req, type_indicator)
    Rlb_ = get_tmp(Rlb, type_indicator)
    Rub_ = get_tmp(Rub, type_indicator)
    Xeq_ = get_tmp(Xeq, type_indicator)
    Xlb_ = get_tmp(Xlb, type_indicator)
    Xub_ = get_tmp(Xub, type_indicator)

    g_update_kernel!(
        PQ,
        a_rhs_,
        v_rhs_,
        Ipeq_,
        Iqeq_,
        Ppf_,
        Qpf_,
        Req_,
        Rlb_,
        Rub_,
        Xeq_,
        Xlb_,
        Xub_,
        dae_t,
        p0_,
        p2i_,
        p2p_,
        p2z_,
        q0_,
        q2i_,
        q2q_,
        q2z_,
        u_,
        v_,
        vcmp_zi_,
        vcmp_zl_,
        vcmp_zu_,
    )
    nothing
end


function g_update_kernel!(
    ::Type{PQ},
    a_rhs,
    v_rhs,
    Ipeq,
    Iqeq,
    Ppf,
    Qpf,
    Req,
    Rlb,
    Rub,
    Xeq,
    Xlb,
    Xub,
    dae_t,
    p0,
    p2i,
    p2p,
    p2z,
    q0,
    q2i,
    q2q,
    q2z,
    u,
    v,
    vcmp_zi,
    vcmp_zl,
    vcmp_zu,
)
    @. a_rhs =
        u * (Ipeq * p2i * v + Ppf * p2p + Req * p2z * v .^ 2) * (dae_t > 0) +
        u * (Rlb * v .^ 2 * vcmp_zl + Rub * v .^ 2 * vcmp_zu + p0 * vcmp_zi) * (dae_t <= 0)
    @. v_rhs =
        u * (Iqeq * q2i * v + Qpf * q2q + Xeq * q2z * v .^ 2) * (dae_t > 0) +
        u * (Xlb * v .^ 2 * vcmp_zl + Xub * v .^ 2 * vcmp_zu + q0 * vcmp_zi) * (dae_t <= 0)
    nothing
end


@views function collect_f!(dae, model::PQ, addr::PQAddr, type_ind)
    nothing
end


@views function collect_g!(dae, model::PQ, addr::PQAddr, type_ind)
    (; a_rhs, v_rhs) = model.rhs_g
    (; a_addr, v_addr) = addr.rhs_g
    get_tmp(dae.g, type_ind)[a_addr] .+= get_tmp(a_rhs, type_ind)
    get_tmp(dae.g, type_ind)[v_addr] .+= get_tmp(v_rhs, type_ind)
    nothing
end


@views function deliver_x!(dae, model::PQ, addr::PQAddr, type_ind)
    nothing
end


@views function deliver_y!(dae, model::PQ, addr::PQAddr, type_ind)
    (; a, v) = model.algeb
    (; a_addr, v_addr) = addr.algeb
    get_tmp(a, type_ind) .= get_tmp(dae.y, type_ind)[a_addr]
    get_tmp(v, type_ind) .= get_tmp(dae.y, type_ind)[v_addr]
    nothing
end


function fx_update!(model::PQ, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function fy_update!(model::PQ, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function gx_update!(model::PQ, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function gy_update!(model::PQ, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    algeb = model.algeb
    config = model.config
    discrete = model.discrete
    jac_gy = model.jac_gy
    param = model.param
    service = model.service

    (; dae_t) = sim_config
    (; v) = algeb
    (; p2i, p2z, q2i, q2z) = config
    (; vcmp_zl, vcmp_zu) = discrete
    (; _gy1, _gy2) = jac_gy
    (; u) = param
    (; Ipeq, Iqeq, Req, Rlb, Rub, Xeq, Xlb, Xub) = service

    gy_update_kernel!(
        PQ,
        _gy1,
        _gy2,
        Ipeq,
        Iqeq,
        Req,
        Rlb,
        Rub,
        Xeq,
        Xlb,
        Xub,
        dae_t,
        p2i,
        p2z,
        q2i,
        q2z,
        u,
        v,
        vcmp_zl,
        vcmp_zu,
    )
    nothing
end


function gy_update_kernel!(
    ::Type{PQ},
    _gy1,
    _gy2,
    Ipeq,
    Iqeq,
    Req,
    Rlb,
    Rub,
    Xeq,
    Xlb,
    Xub,
    dae_t,
    p2i,
    p2z,
    q2i,
    q2z,
    u,
    v,
    vcmp_zl,
    vcmp_zu,
)
    @. _gy1 =
        u * (Ipeq * p2i + 2 * Req * p2z * v) * (dae_t > 0) +
        u * (2 * Rlb * v * vcmp_zl + 2 * Rub * v * vcmp_zu) * (dae_t <= 0)
    @. _gy2 =
        u * (Iqeq * q2i + 2 * Xeq * q2z * v) * (dae_t > 0) +
        u * (2 * Xlb * v * vcmp_zl + 2 * Xub * v * vcmp_zu) * (dae_t <= 0)
    nothing
end


function l_update_var!(model::PQ; dae_t = 0.0)
    # --- expand structs ---
    algeb = model.algeb
    discrete = model.discrete
    param = model.param

    (; v) = algeb
    (; vcmp_zi, vcmp_zl, vcmp_zu) = discrete
    (; vmin, vmax) = param

    # --- begin `l_update_var` calls ---
    l_update_var!(
        Val{:Limiter},
        v,
        vmin,
        vmax,
        vcmp_zi,
        vcmp_zl,
        vcmp_zu;
        equal = false,
        enable = true,
        no_lower = false,
        no_upper = false,
        sign_lower = 1,
        sign_upper = 1,
    )
    nothing
end


function l_update_eq!(model::PQ; dae_t = 0.0)
    # empty function
    nothing
end
