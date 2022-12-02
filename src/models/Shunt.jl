"""
Shunt Union-type parameters
"""
struct ShuntUnionParam
    idx::Vector{String}
    name::Vector{String}
    bus::Vector{String}
end

"""
Shunt parameter definitions
"""
struct ShuntParam{T}
    u::T
    Sn::T
    Vn::T
    g::T
    b::T
    fn::T
end

struct ShuntState{T} end


ShuntState() = ShuntState{Nothing}()


struct ShuntStateAddr end


struct ShuntAlgeb{T}
    a::T
    v::T
end


struct ShuntAlgebAddr{T<:AbstractArray}
    a_addr::T
    v_addr::T
end


struct ShuntRhsF{T} end


ShuntRhsF() = ShuntRhsF{Nothing}()


struct ShuntRhsFAddr end


struct ShuntRhsG{T}
    a_rhs::T
    v_rhs::T
end


struct ShuntRhsGAddr{T<:AbstractArray}
    a_addr::T
    v_addr::T
end


struct ShuntService{} end


struct ShuntDiscrete{T} end


ShuntDiscrete() = ShuntDiscrete{Nothing}()


struct ShuntConfig
    allow_adjust::Float64
    adjust_lower::Float64
    adjust_upper::Float64
end


struct ShuntJacFx{T} end


struct ShuntJacFy{T} end


struct ShuntJacGx{T} end


struct ShuntJacGy{T}
    _gy1::T
    _gy2::T
end


"""

    Phasor-domain shunt compensator Model.
    
"""
struct Shunt{P,S,A,F,G,SV,D}
    n::Int64
    uparam::ShuntUnionParam
    param::ShuntParam{P}
    state::ShuntState{S}
    algeb::ShuntAlgeb{A}
    rhs_f::ShuntRhsF{F}
    rhs_g::ShuntRhsG{G}
    service::SV
    discrete::D
    config::ShuntConfig
end


struct ShuntAddr
    state::ShuntStateAddr
    algeb::ShuntAlgebAddr
    rhs_f::ShuntRhsFAddr
    rhs_g::ShuntRhsGAddr
end


function f_update!(
    model::Shunt,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function g_update!(
    model::Shunt,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    param = model.param
    rhs_g = model.rhs_g

    (; v) = algeb
    (; b, g, u) = param
    (; a_rhs, v_rhs) = rhs_g
    v_ = get_tmp(v, type_indicator)
    b_ = get_tmp(b, type_indicator)
    g_ = get_tmp(g, type_indicator)
    u_ = get_tmp(u, type_indicator)
    a_rhs_ = get_tmp(a_rhs, type_indicator)
    v_rhs_ = get_tmp(v_rhs, type_indicator)

    g_update_kernel!(Shunt, a_rhs_, v_rhs_, b_, g_, u_, v_)
    nothing
end


function g_update_kernel!(::Type{Shunt}, a_rhs, v_rhs, b, g, u, v)
    @. a_rhs = g * u * v .^ 2
    @. v_rhs = -b * u * v .^ 2
    nothing
end


@views function collect_f!(dae, model::Shunt, addr::ShuntAddr, type_ind)
    nothing
end


@views function collect_g!(dae, model::Shunt, addr::ShuntAddr, type_ind)
    (; a_rhs, v_rhs) = model.rhs_g
    (; a_addr, v_addr) = addr.rhs_g
    get_tmp(dae.g, type_ind)[a_addr] .+= get_tmp(a_rhs, type_ind)
    get_tmp(dae.g, type_ind)[v_addr] .+= get_tmp(v_rhs, type_ind)
    nothing
end


@views function deliver_x!(dae, model::Shunt, addr::ShuntAddr, type_ind)
    nothing
end


@views function deliver_y!(dae, model::Shunt, addr::ShuntAddr, type_ind)
    (; a, v) = model.algeb
    (; a_addr, v_addr) = addr.algeb
    get_tmp(a, type_ind) .= get_tmp(dae.y, type_ind)[a_addr]
    get_tmp(v, type_ind) .= get_tmp(dae.y, type_ind)[v_addr]
    nothing
end


function fx_update!(
    model::Shunt,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function fy_update!(
    model::Shunt,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gx_update!(
    model::Shunt,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gy_update!(
    model::Shunt,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    algeb = model.algeb
    jac_gy = model.jac_gy
    param = model.param

    (; v) = algeb
    (; _gy1, _gy2) = jac_gy
    (; b, g, u) = param

    gy_update_kernel!(Shunt, _gy1, _gy2, b, g, u, v)
    nothing
end


function gy_update_kernel!(::Type{Shunt}, _gy1, _gy2, b, g, u, v)
    @. _gy1 = 2 * g * u * v
    @. _gy2 = -2 * b * u * v
    nothing
end


function l_update_var!(model::Shunt; dae_t = 0.0)
    # empty function
    nothing
end


function l_update_eq!(model::Shunt; dae_t = 0.0)
    # empty function
    nothing
end
