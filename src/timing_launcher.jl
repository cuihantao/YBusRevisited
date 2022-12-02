include("scripts/base_import.jl");


if @isdefined(MATPOWER_DATA_PATH) == false
    error(
        "`MATPOWER_DATA_PATH` not defined. Define this variable to point to the matpower data folder.",
    )
end

BLAS.get_config()
versioninfo()

casenames = [
    "case14.m",
    "case118.m",
    "case300.m",
    "case1354pegase.m",
    "case2736sp.m",
    "case9241pegase.m",
    "case_ACTIVSg25k.m",
    "case_ACTIVSg70k.m",
    "case_SyntheticUSA.m",
];

# define a global variable for the case
global case = "";


bresults_ewise_cases = Dict{String,Vector{BenchmarkTools.Trial}}();
bresults_Ybus_cases = Dict{String,Vector{BenchmarkTools.Trial}}();

for casename in casenames
    global case = joinpath(MATPOWER_DATA_PATH, casename)

    # ------------------------------------------------
    # Benchmark element-wise method
    # ------------------------------------------------

    include("scripts/timing_elemwise.jl")

    # ------------------------------------------------
    # Benchmark Ybus method
    # ------------------------------------------------
    include("scripts/timing_Ybus.jl")

    # ------------------------------------------------
    # Cross-verify results
    # ------------------------------------------------

    include("scripts/verify_correctness.jl")

    bresults_Ybus_cases[casename] = bresults_Ybus
    bresults_ewise_cases[casename] = bresults_ewise

end


@show bresults_ewise_cases

@show bresults_Ybus_cases

#=
# Uncomment block for serialization

using Serialization
using Base.Sys

serialize(joinpath("benchmark_data/", String(Base.Sys.cpu_info()[1].model)),
          (bresults_ewise_cases, bresults_Ybus_cases))
=#
