module Snowflakes

using XLSX
using DataFrames
using JuMP
using Query
using Distances
using Random
using UUIDs
using MathOptInterface
using Combinatorics
using GLPK

include("Structures.jl")
include("Base.jl")

export base
export stats
export initStab
export root

end
