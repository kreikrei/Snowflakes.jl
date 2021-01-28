using Snowflakes
using Test
using JuMP
using GLPK

path = joinpath(@__DIR__,"testdata2.xlsx")
extract!(path)

base()
using Gurobi
set_default_optimizer!(Gurobi.Optimizer)

test_root = root(;slC=1000.0,suC=-1000.0)

test_mp = master(test_root;silent=false)
test_dual = Snowflakes.getDuals(test_mp)

test_sp = Snowflakes.sub(test_root,test_dual;silent=false)

test_x = value.(test_sp.obj_dict[:x])
test_p = value.(test_sp.obj_dict[:p])

test_x.data
matrix_x(k,f,t) = JuMP.Containers.DenseAxisArray{Float64}(undef,base().K[k].cover,base().K[k].cover)
matrix_x(1,1,1)


getTours(test_x)


test_col = Snowflakes.getCols(test_sp)
push!(test_root.columns,test_col)

println(value.(test_mp.obj_dict[:θ]))




@testset "Base.jl" begin
    @test stats().number_of_vertices == 46
    @test stats().number_of_vehicles == 84

    @test isequal( #all in cover_list is in keys(V)
        sort(stats().cover_list),sort(collect(keys(base().V)))
    )

    idx = rand(collect(keys(base().V))) #random point
    @test base().dist[idx,idx] == 999999999

    @test initStab(;slC = 500.0, suC = -500.0).slLim == abs.(base().d)
    @test initStab(;slC = 500.0, suC = -500.0).suLim == abs.(base().d)

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

    @test has_values(test_mp)
end

using Gurobi
set_default_optimizer!(Gurobi.Optimizer)
test_root = root(;slC=500.0,suC=-500.0)
test_mp = master(test_root;silent=false)
test_duals = Snowflakes.getDuals(test_mp)
test_sub = Snowflakes.sub(test_root,test_duals;silent=false)
