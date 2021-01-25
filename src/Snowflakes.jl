module Snowflakes

#core opt
using JuMP
using MathOptInterface

#data extraction
using XLSX
using DataFrames
using Query
using Distances

#support
using Random
using Combinatorics
using UUIDs

include("Structures.jl")
include("Base.jl")
include("Column.jl")

export base
export stats
export initStab
export root

end
