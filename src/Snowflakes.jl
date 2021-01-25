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
include("Settings.jl")

export extract
export base!
export stats
export initStab
export root

export set_default_optimizer!
export get_default_optimizer
export reset_default_optimizer

export master

end
