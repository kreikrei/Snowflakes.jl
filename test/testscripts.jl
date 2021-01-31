using Snowflakes
using Test
using JuMP
using GLPK

#INITIALIZE DATA
path = joinpath(@__DIR__,"testdata2.xlsx")
extract!(path)

@testset "Base.jl" begin
    @test isequal( #all in cover_list is in keys(V)
        sort(stats().cover_list),sort(collect(keys(b().V)))
    )

    idx = rand(collect(keys(b().V))) #random point d[i,i]
    @test b().dist[idx,idx] == 999999999

    @test initStab(;slC = 500.0, suC = -500.0).slLim == abs.(b().d)
    @test initStab(;slC = 500.0, suC = -500.0).suLim == abs.(b().d)

    @test isa(root(;slC = 500.0, suC = -500.0),Snowflakes.node)
end

@testset "Settings.jl" begin
    set_default_optimizer!(GLPK.Optimizer)

    @test get_default_optimizer() == GLPK.Optimizer
    @test reset_default_optimizer() == nothing
end

@testset "Column.jl" begin
    set_default_optimizer!(GLPK.Optimizer)
    test_root = root(;slC=500.0,suC=-500.0)
    test_mp = master(test_root;silent=false)

    @test has_values(test_mp) #has primal values
    @test has_duals(test_mp) #has dual values

    test_duals = getDuals(test_mp)
    test_sp = sub(test_root,test_duals;silent=false)

    @test objective_value(test_sp) < 0 #first iter sub < 0
    @test updateStab!(test_root.stab,0.5).slLim == 0.5 .* b().d
end
