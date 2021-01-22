using Snowflakes
using Test

@testset "Snowflakes.jl" begin
    path = joinpath(@__DIR__,"testdata.xlsx")
    res = base(path)

    @test report(res).number_of_vertices == 46
    @test report(res).number_of_vehicles == 84

    @test isequal(
        sort(report(res).unique_cover),sort(collect(keys(res.V)))
    )

    @test res.dist[1,1] == 999999999
end
