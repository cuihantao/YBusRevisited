include("struct_map.jl");

using PreallocationTools

Indicator(x) = x;


"""
Create an instance of an YBusRevisited model from an `andes.System`.
"""
function load_model(ss, model_name::Symbol)

    jl_struct = model_map[model_name]
    model_name_str = string(model_name)
    ss_model = ss.__dict__[model_name_str]
    inputs = ss_model._input

    # === create empty instances for models with no device ===
    if ss_model.n == 0
        return jl_struct()
    end

    members = Vector{Any}()

    push!(members, ss_model.n)

    for ftype in fields_map[model_name]

        fvalues = Vector{Any}()
        type_str = string(ftype)

        if occursin("Param", type_str) || occursin("Service", type_str)

            fnames = fieldnames(ftype)

            for name in fnames

                val = zeros(ss_model.n)

                if occursin("_C", string(name))
                    val = zeros(ComplexF64, ss_model.n)
                end

                if hasproperty(ss_model, name)
                    val = getproperty(ss_model, name).v
                end

                if occursin("UnionParam", type_str)
                    val = string.(val)
                end

                push!(fvalues, val)
            end

            soa_view = ftype(fvalues...)
            push!(members, soa_view)

        elseif occursin("State", type_str) || occursin("Algeb", type_str)
            fnames = fieldnames(ftype)

            for name in fnames
                val = getproperty(ss_model, name).v
                val = dualcache(val)
                push!(fvalues, val)
            end

            soa_view = ftype(fvalues...)
            push!(members, soa_view)

        elseif occursin("Rhs", type_str)
            fnames = fieldnames(ftype)
            for name in fnames
                varname = replace(String(name), "_rhs" => "")
                val = getproperty(ss_model, varname).e
                val = dualcache(val)

                push!(fvalues, val)
            end

            soa_view = ftype(fvalues...)
            push!(members, soa_view)

            # elseif occursin("Jac", type_str)

            #     fnames = fieldnames(ftype)

            #     for name in fnames
            #         push!(fvalues, [0])
            #     end

            #     soa_view = ftype(fvalues...)
            #     push!(members, soa_view)

        elseif occursin("Discrete", type_str)
            fnames = fieldnames(ftype)
            for name in fnames
                val = inputs[string(name)]
                push!(fvalues, val)
            end

            soa_view = ftype(fvalues...)
            push!(members, soa_view)

        elseif occursin("Config", type_str)
            fnames = fieldnames(ftype)
            for name in fnames
                push!(fvalues, getproperty(ss_model.config, name))
            end
            push!(members, ftype(fvalues...))
        end
    end

    return jl_struct(members...)
end


"""
Create a system using a dictionary of {Symbol, Model}.

Parameters
----------
skip_empty : bool, default to true
    True to skip the creation of models with zero devices

Returns
-------
Dict : Symbol => Model
    Model name symbol mapped to instances

"""
function load_sys_data(ss, skip_empty = true)

    sys_data = Dict{Symbol,Any}()
    for fname in keys(model_map)
        if (skip_empty == true) && (getproperty(ss, fname).n == 0)
            continue
        end

        sys_data[fname] = load_model(ss, fname)
    end

    return sys_data
end


"""
Load an ANDES case into a Julia Dict.

Returns a Julia Dict and a PyObject for andes.System.
"""
function load_system(ss)
    sys_data = load_sys_data(ss)
    models = collect(keys(sys_data))
    sys_config = load_sys_config(ss)
    sim_config = load_sim_config(ss)
    dae = load_dae(ss)
    addr = load_address(ss)

    return System(models, sys_data, sys_config, dae, addr, sim_config)
end


"""
Helper function to load variable and equation addresses for a model.
"""
function _load_address_for_model(ss, model_name)

    out = Vector{Any}()

    model_name_str = string(model_name)
    ss_model = ss.__dict__[model_name_str]

    addr_types = address_substructs[model_name]

    for addr_type in addr_types

        fvalues = Vector{Vector}()

        fnames = fieldnames(addr_type)

        for name in fnames
            name = replace(String(name), "_addr" => "")

            push!(fvalues, getproperty(ss_model, name).a .+ 1)  # convert to 1-based indices
        end

        push!(out, addr_type(fvalues...))
    end

    return address_map[model_name](out...)

end


"""
Load variable and equation addresses from an ANDES System.
"""
function load_address(ss, skip_empty = true)

    addresses = Dict{Symbol,Any}()

    for modelName in keys(model_map)
        if (skip_empty == true) && (getproperty(ss, modelName).n == 0)
            continue
        end

        addresses[modelName] = _load_address_for_model(ss, modelName)
    end

    return addresses
end


"""
Load numerical DAE structure.
"""
function load_dae(ss)

    fvalues = Vector()
    for fname in fieldnames(DAE)

        val = getproperty(ss.dae, fname)

        # use PreallocationTools.DiffCache to support autodiff
        if isa(val, AbstractArray)
            val = dualcache(val)
        end

        push!(fvalues, val)
    end

    return DAE(fvalues...)
end


"""
Load config for system

TODO: implement data loading from system.
"""
function load_sys_config(ss)

    return SysConfig(100.0, 60.0)
end


function load_sim_config(ss)

    return SimConfig(0.0)
end
