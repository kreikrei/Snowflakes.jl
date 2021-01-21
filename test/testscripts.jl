using Snowflakes
using Test

@testset "Snowflakes.jl" begin
    @test report(base(joinpath(@__DIR__,"testdata.xlsx"))).number_of_vertices == 46
    @test report(base(joinpath(@__DIR__,"testdata.xlsx"))).number_of_vehicles == 84
end
