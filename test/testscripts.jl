using Snowflakes
using Test
using JuMP
using GLPK
using Gurobi
using DataFrames

#SET SOLVER
GUROBI_ENV = Gurobi.Env()
set_default_optimizer!(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))

#INITIALIZE DATA
path = joinpath(@__DIR__,"testdata1.xlsx")
extract!(path)

#GENERATE ROOT
test_root = root(;slC=1000.0,suC=-1000.0)

#colGen
colGen(test_root;silent=true,track=true,maxCG=Inf)

test_z = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
test_R = Dict(1:length(test_root.columns) .=> test_root.columns)
test_mp = master(test_root;silent=false)
θ = test_mp.obj_dict[:θ]
collection = DataFrame(i=Int64[],k=Int64[],t=Int64[],val=Float64[])

for i in keys(b().V), k in keys(b().K), t in b().T
    iter = keys(filter(p -> last(p).z[i,k,t] > 0,test_R))

    if !isempty(iter)
        tot = value(sum(
            θ[r,k,t]
            for r in iter
        ))

        if tot > 0 && !isinteger(tot)
            append!(collection,DataFrame(i=i,k=k,t=t,val=tot))
        end
    end
end

push!(test_root.bounds,Snowflakes.bound((i=15,k=72,t=3),1))

master(test_root;silent=false)

println(collection)

test_z[25,43,1]

colGen(test_root;silent=true,maxCG=10.0,track=true)





@testset "Base.jl" begin
    @test stats().number_of_vertices == 46
    @test stats().number_of_vehicles == 84

    @test isequal( #all in cover_list is in keys(V)
        sort(stats().cover_list),sort(collect(keys(b().V)))
    )

    idx = rand(collect(keys(b().V))) #random point
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

    @test has_values(test_mp)
end
