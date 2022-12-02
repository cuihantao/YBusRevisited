using YBusRevisited
using Test

@testset "YBusRevisited.jl" begin
    include(joinpath(@__DIR__, "test_load.jl"))
    include(joinpath(@__DIR__, "test_residual_match.jl"))
end
