import PreallocationTools: get_tmp

get_tmp(dc::AbstractArray, u) = dc;

get_tmp(dc::Float64, u) = dc;


struct SysConfig
    sys_mva::Float64
    sys_f::Float64
end


struct SimConfig
    dae_t::Float64
end


struct System
    models::Vector{Symbol}
    data::Dict{Symbol,Any}
    config::SysConfig
    dae::DAE
    addr::Dict{Symbol,Any}
    sim_config::SimConfig
end


"""
$(SIGNATURES)

API function to update selected or all residual equations.
"""
function update_residuals!(sys::System, type_indicator, eqn_name::Symbol = :fg)

    eqn_name in (:f, :g, :fg) || error("Unexpected equation name $eqn_name")

    if eqn_name in (:f, :fg)
        YBusRevisited.f_update!(sys, type_indicator)
    end

    if eqn_name in (:g, :fg)
        YBusRevisited.g_update!(sys, type_indicator)
    end

    nothing
end


"""
$(SIGNATURES)

API function to collect all algebraic residual values to `dae.g`.
"""
function collect_residuals!(sys::System, type_indicator, eqn_name::Symbol = :fg)
    eqn_name in (:f, :g, :fg) || error("Unexpected equation name $eqn_name")

    if eqn_name in (:f, :fg)
        YBusRevisited.collect_f!(sys, type_indicator)
    end

    if eqn_name in (:g, :fg)
        YBusRevisited.collect_g!(sys, type_indicator)
    end

    nothing
end


"""
API function for delivering variables to model data.
"""
function deliver_variables!(sys::System, type_indicator, var_name::Symbol = :xy)
    var_name in (:x, :y, :xy) || error("Unexpected variable name $var_name")

    if var_name in (:x, :xy)
        YBusRevisited.deliver_x!(sys, type_indicator)
    end

    if var_name in (:x, :xy)
        YBusRevisited.deliver_y!(sys, type_indicator)
    end

    nothing
end


function deliver_x!(sys::System, type_indicator)
    for model_name in sys.models
        YBusRevisited.deliver_x!(
            sys.dae,
            sys.data[model_name],
            sys.addr[model_name],
            type_indicator,
        )
    end

    nothing

end


function deliver_y!(sys::System, type_indicator)
    for model_name in sys.models
        YBusRevisited.deliver_y!(
            sys.dae,
            sys.data[model_name],
            sys.addr[model_name],
            type_indicator,
        )
    end

    nothing
end


"""
Internal wrapper function to evaluate residual equations and store results to
model internal field `rhs_g`.
"""
function g_update!(sys::System, type_indicator)

    for model_name in sys.models
        YBusRevisited.g_update!(
            sys.data[model_name],
            sys.config,
            sys.sim_config,
            type_indicator,
        )
    end

    nothing
end


"""
Internal wrapper function to evaluate differential RHS equations and store
results to model internal field `rhs_f`.
"""
function f_update!(sys::System, type_indicator)

    for model_name in sys.models
        YBusRevisited.f_update!(
            sys.data[model_name],
            sys.config,
            sys.sim_config,
            type_indicator,
        )
    end

    nothing
end



"""
Internal wrapper function to collect all algebraic equation residual values to
`dae.g`.
"""
function collect_g!(sys::System, type_indicator)

    for model_name in sys.models
        YBusRevisited.collect_g!(
            sys.dae,
            sys.data[model_name],
            sys.addr[model_name],
            type_indicator,
        )
    end

    nothing
end


"""
Internal wrapper function to collect all differential equation residual values
to `dae.f`.
"""
function collect_f!(sys::System, type_indicator)

    for model_name in sys.models
        YBusRevisited.collect_f!(
            sys.dae,
            sys.data[model_name],
            sys.addr[model_name],
            type_indicator,
        )
    end

    nothing
end


"""
Evaluate and in-place update the system residual array for the given variable
inputs.
"""
@views function evaluate_residuals!(
    sys::YBusRevisited.System,
    out::AbstractArray,
    xy::AbstractArray,
)

    dae_x = get_tmp(sys.dae.x, xy)
    dae_y = get_tmp(sys.dae.y, xy)

    dae_x .= xy[1:sys.dae.n]
    dae_y .= xy[sys.dae.n+1:sys.dae.n+sys.dae.m]

    YBusRevisited.deliver_variables!(sys, xy)

    get_tmp(sys.dae.f, xy) .= 0
    get_tmp(sys.dae.g, xy) .= 0

    YBusRevisited.update_residuals!(sys, xy)
    YBusRevisited.collect_residuals!(sys, xy)

    @views out[1:sys.dae.n] .= get_tmp(sys.dae.f, xy)
    @views out[sys.dae.n+1:sys.dae.m+sys.dae.n] = get_tmp(sys.dae.g, xy)

    nothing
end


"""
Helper function to return a closure that evaluates the system residuals and
updates the output array in place, using the input variable values.

The closure uses the signature `f!(output, input)`, which is compatible with
DiffEq solvers.
"""
function get_residual_function(sys::YBusRevisited.System)

    return function res!(out::AbstractArray, xy::AbstractArray)
        evaluate_residuals!(sys, out, xy)
    end

end
