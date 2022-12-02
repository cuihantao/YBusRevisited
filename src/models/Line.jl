"""
Line Union-type parameters
"""
struct LineUnionParam
    idx::Vector{String}
    name::Vector{String}
    bus1::Vector{String}
    bus2::Vector{String}
    owner::Vector{String}
    xcoord::Vector{String}
    ycoord::Vector{String}
end

"""
Line parameter definitions
"""
struct LineParam{T}
    u::T
    Sn::T
    fn::T
    Vn1::T
    Vn2::T
    r::T
    x::T
    b::T
    g::T
    b1::T
    g1::T
    b2::T
    g2::T
    trans::T
    tap::T
    phi::T
    rate_a::T
    rate_b::T
    rate_c::T
end

struct LineState{T} end


LineState() = LineState{Nothing}()


struct LineStateAddr end


struct LineAlgeb{T}
    a1::T
    a2::T
    v1::T
    v2::T
end


struct LineAlgebAddr{T<:AbstractArray}
    a1_addr::T
    a2_addr::T
    v1_addr::T
    v2_addr::T
end


struct LineRhsF{T} end


LineRhsF() = LineRhsF{Nothing}()


struct LineRhsFAddr end


struct LineRhsG{T}
    a1_rhs::T
    a2_rhs::T
    v1_rhs::T
    v2_rhs::T
end


struct LineRhsGAddr{T<:AbstractArray}
    a1_addr::T
    a2_addr::T
    v1_addr::T
    v2_addr::T
end


struct LineService{C,T}
    gh::T
    bh::T
    gk::T
    bk::T
    yh::C
    yk::C
    yhk::C
    ghk::T
    bhk::T
    itap::T
    itap2::T
    bhbhk::T
    ghghk::T
    sine::T
    cosine::T
    v1v2::T
    v1_2::T
    v2_2::T
    bhksine::T
    bhkcosine::T
    ghksine::T
    ghkcosine::T
    itapv1v2::T
    itapv1::T
    itapv2::T
end


struct LineDiscrete{T} end


LineDiscrete() = LineDiscrete{Nothing}()


struct LineConfig
    allow_adjust::Float64
    adjust_lower::Float64
    adjust_upper::Float64
end




struct LineJacFx{T} end


struct LineJacFy{T} end


struct LineJacGx{T} end


struct LineJacGy{T}
    _gy1::T
    _gy2::T
    _gy3::T
    _gy4::T
    _gy5::T
    _gy6::T
    _gy7::T
    _gy8::T
    _gy9::T
    _gy10::T
    _gy11::T
    _gy12::T
    _gy13::T
    _gy14::T
    _gy15::T
    _gy16::T
end


"""

    AC transmission line model.

    The model is also used for two-winding transformer. Transformers can set the
    tap ratio in ``tap`` and/or phase shift angle ``phi``.

    To reduce the number of variables, line injections are summed at bus
    equations and are not stored. Current injections are not computed.

"""
struct Line{P,S,A,F,G,SV,D}
    n::Int64
    uparam::LineUnionParam
    param::LineParam{P}
    state::LineState{S}
    algeb::LineAlgeb{A}
    rhs_f::LineRhsF{F}
    rhs_g::LineRhsG{G}
    service::SV
    discrete::D
    config::LineConfig
end


struct LineAddr
    state::LineStateAddr
    algeb::LineAlgebAddr
    rhs_f::LineRhsFAddr
    rhs_g::LineRhsGAddr
end


function f_update!(
    model::Line,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function g_update!(
    model::Line,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    param = model.param
    rhs_g = model.rhs_g
    service = model.service

    (; a1, a2, v1, v2) = algeb
    (; phi, u) = param
    (; a1_rhs, a2_rhs, v1_rhs, v2_rhs) = rhs_g
    (; bh, bhk, gh, ghk, itap, itap2) = service
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

    g_update_kernel!(
        Line,
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
    )
    nothing
end


function g_update_kernel!(
    ::Type{Line},
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
)
    @. a1_rhs =
        u * (
            -itap * v1 * v2 * (-bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi)) +
            itap2 * v1 .^ 2 * (gh + ghk)
        )
    @. a2_rhs =
        u * (
            -itap * v1 * v2 * (bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi)) +
            v2 .^ 2 * (gh + ghk)
        )
    @. v1_rhs =
        u * (
            -itap * v1 * v2 * (-bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi)) -
            itap2 * v1 .^ 2 * (bh + bhk)
        )
    @. v2_rhs =
        u * (
            itap * v1 * v2 * (bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi)) -
            v2 .^ 2 * (bh + bhk)
        )
    nothing
end


@views function collect_f!(dae, model::Line, addr::LineAddr, type_ind)
    nothing
end


@views function collect_g!(dae, model::Line, addr::LineAddr, type_ind)
    (; a1_rhs, a2_rhs, v1_rhs, v2_rhs) = model.rhs_g
    (; a1_addr, a2_addr, v1_addr, v2_addr) = addr.rhs_g
    get_tmp(dae.g, type_ind)[a1_addr] .+= get_tmp(a1_rhs, type_ind)
    get_tmp(dae.g, type_ind)[a2_addr] .+= get_tmp(a2_rhs, type_ind)
    get_tmp(dae.g, type_ind)[v1_addr] .+= get_tmp(v1_rhs, type_ind)
    get_tmp(dae.g, type_ind)[v2_addr] .+= get_tmp(v2_rhs, type_ind)
    nothing
end


@views function deliver_x!(dae, model::Line, addr::LineAddr, type_ind)
    nothing
end


@views function deliver_y!(dae, model::Line, addr::LineAddr, type_ind)
    (; a1, a2, v1, v2) = model.algeb
    (; a1_addr, a2_addr, v1_addr, v2_addr) = addr.algeb
    get_tmp(a1, type_ind) .= get_tmp(dae.y, type_ind)[a1_addr]
    get_tmp(a2, type_ind) .= get_tmp(dae.y, type_ind)[a2_addr]
    get_tmp(v1, type_ind) .= get_tmp(dae.y, type_ind)[v1_addr]
    get_tmp(v2, type_ind) .= get_tmp(dae.y, type_ind)[v2_addr]
    nothing
end


function fx_update!(
    model::Line,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function fy_update!(
    model::Line,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gx_update!(
    model::Line,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gy_update!(
    model::Line,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    jac_gy = model.jac_gy
    param = model.param
    service = model.service

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
    (; bh, bhk, gh, ghk, itap, itap2) = service

    gy_update_kernel!(
        Line,
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
        v1,
        v2,
    )
    nothing
end



function gy_update_kernel!(
    ::Type{Line},
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
    v1,
    v2,
)
    @. _gy1 = -itap * u * v1 * v2 * (bhk * cos(-a1 + a2 + phi) + ghk * sin(-a1 + a2 + phi))
    @. _gy2 = -itap * u * v1 * v2 * (-bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi))
    @. _gy3 =
        u * (
            -itap * v2 * (-bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi)) +
            2 * itap2 * v1 * (gh + ghk)
        )
    @. _gy4 = -itap * u * v1 * (-bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi))
    @. _gy5 = -itap * u * v1 * v2 * (-bhk * cos(-a1 + a2 + phi) + ghk * sin(-a1 + a2 + phi))
    @. _gy6 = -itap * u * v1 * v2 * (bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi))
    @. _gy7 = -itap * u * v2 * (bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi))
    @. _gy8 =
        u * (
            -itap * v1 * (bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi)) +
            2 * v2 * (gh + ghk)
        )
    @. _gy9 = -itap * u * v1 * v2 * (-bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi))
    @. _gy10 = -itap * u * v1 * v2 * (bhk * sin(-a1 + a2 + phi) - ghk * cos(-a1 + a2 + phi))
    @. _gy11 =
        u * (
            -itap * v2 * (-bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi)) -
            2 * itap2 * v1 * (bh + bhk)
        )
    @. _gy12 = -itap * u * v1 * (-bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi))
    @. _gy13 = itap * u * v1 * v2 * (bhk * sin(-a1 + a2 + phi) + ghk * cos(-a1 + a2 + phi))
    @. _gy14 = itap * u * v1 * v2 * (-bhk * sin(-a1 + a2 + phi) - ghk * cos(-a1 + a2 + phi))
    @. _gy15 = itap * u * v2 * (bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi))
    @. _gy16 =
        u * (
            itap * v1 * (bhk * cos(-a1 + a2 + phi) - ghk * sin(-a1 + a2 + phi)) -
            2 * v2 * (bh + bhk)
        )
    nothing
end


function l_update_var!(model::Line; dae_t = 0.0)
    # empty function
    nothing
end


function l_update_eq!(model::Line; dae_t = 0.0)
    # empty function
    nothing
end
