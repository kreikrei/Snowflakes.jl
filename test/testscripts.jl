using Snowflakes
using Test
using GLPK

path = joinpath(@__DIR__,"testdata.xlsx")
base!(path)

@testset "Base.jl" begin
    @test stats().number_of_vertices == 46
    @test stats().number_of_vehicles == 84

    @test isequal( #all in cover_list is in keys(V)
        sort(stats().cover_list),sort(collect(keys(extract().V)))
    )

    idx = rand(collect(keys(extract().V))) #random point
    @test extract().dist[idx,idx] == 999999999

    @test initStab(;slC = 500.0, suC = -500.0).slLim == abs.(extract().d)
    @test initStab(;slC = 500.0, suC = -500.0).suLim == abs.(extract().d)

    @test isa(root(;slC = 500.0, suC = -500.0),Snowflakes.node)
end

@testset "Settings.jl" begin
    set_default_optimizer!(GLPK.Optimizer)

    @test get_default_optimizer() == GLPK.Optimizer
    @test reset_default_optimizer() == nothing
end
