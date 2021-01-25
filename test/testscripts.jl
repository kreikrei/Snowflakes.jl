using Snowflakes
using Test

@testset "Base.jl test" begin
    path = joinpath(@__DIR__,"testdata.xlsx")
    res = base(path)

    @test stats(res).number_of_vertices == 46
    @test stats(res).number_of_vehicles == 84

    @test isequal( #all in cover_list is in keys(V)
        sort(stats(res).cover_list),sort(collect(keys(res.V)))
    )

    idx = rand(collect(keys(res.V))) #random point
    @test res.dist[idx,idx] == 999999999

    @test initStab(res;slC = 500.0, suC = -500.0).slLim == abs.(res.d)
    @test initStab(res;slC = 500.0, suC = -500.0).suLim == abs.(res.d)

    @test isa(root(res;slC = 500.0, suC = -500.0),Snowflakes.node)
end
