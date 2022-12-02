"""
Bus Union-type parameters
"""
struct BusUnionParam
    idx::Vector{String}
    name::Vector{String}
    xcoord::Vector{String}
    ycoord::Vector{String}
    area::Vector{String}
    zone::Vector{String}
    owner::Vector{String}
end

"""
Bus parameter definitions
"""
struct BusParam{T}
    u::T
    Vn::T
    vmax::T
    vmin::T
    v0::T
    a0::T
end

struct BusState{T} end


BusState() = BusState{Nothing}()


struct BusStateAddr end


struct BusAlgeb{T}
    a::T
    v::T
end


struct BusAlgebAddr{T<:AbstractArray}
    a_addr::T
    v_addr::T
end


struct BusRhsF{T} end


BusRhsF() = BusRhsF{Nothing}()


struct BusRhsFAddr end


struct BusRhsG{T} end


BusRhsG() = BusRhsG{Nothing}()


struct BusRhsGAddr end


struct BusService{} end


struct BusDiscrete{T} end


BusDiscrete() = BusDiscrete{Nothing}()


struct BusConfig
    allow_adjust::Float64
    adjust_lower::Float64
    adjust_upper::Float64
    flat_start::Float64
end


struct BusJacFx{T} end


struct BusJacFy{T} end


struct BusJacGx{T} end


struct BusJacGy{T} end


"""

    AC Bus model.

    Power balance equation have the form of ``load - injection = 0``.
    Namely, load is positively summed, while injections are negative.
    
"""
struct Bus{P,S,A,F,G,SV,D}
    n::Int64
    uparam::BusUnionParam
    param::BusParam{P}
    state::BusState{S}
    algeb::BusAlgeb{A}
    rhs_f::BusRhsF{F}
    rhs_g::BusRhsG{G}
    service::SV
    discrete::D
    config::BusConfig
end


struct BusAddr
    state::BusStateAddr
    algeb::BusAlgebAddr
    rhs_f::BusRhsFAddr
    rhs_g::BusRhsGAddr
end


function f_update!(model::Bus, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


function g_update!(model::Bus, sys_config::SysConfig, sim_config::SimConfig, type_indicator)
    nothing
end


@views function collect_f!(dae, model::Bus, addr::BusAddr, type_ind)
    nothing
end


@views function collect_g!(dae, model::Bus, addr::BusAddr, type_ind)
    nothing
end


@views function deliver_x!(dae, model::Bus, addr::BusAddr, type_ind)
    nothing
end


@views function deliver_y!(dae, model::Bus, addr::BusAddr, type_ind)
    nothing
end


function fx_update!(
    model::Bus,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function fy_update!(
    model::Bus,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gx_update!(
    model::Bus,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function gy_update!(
    model::Bus,
    sys_config::SysConfig,
    sim_config::SimConfig,
    type_indicator,
)
    nothing
end


function l_update_var!(model::Bus; dae_t = 0.0)
    # empty function
    nothing
end


function l_update_eq!(model::Bus; dae_t = 0.0)
    # empty function
    nothing
end
