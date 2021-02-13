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
using Statistics

include("Structures.jl")
include("Base.jl")
include("Column.jl")
include("Settings.jl")
include("Branching.jl")

export set_default_optimizer!
export get_default_optimizer
export reset_default_optimizer
export set_slack_coeff!
export sl_C
export set_surp_coeff!
export su_C
export set_silence!
export silent

export b
export qmax
export imax
export extract!
export stats
export initStab
export root
export origin

export Q
export f
export tot
export master
export buildMaster
export getDuals
export sub
export buildSub
export getCols
export updateStab!
export checkStab
export colGen

export separate
export vtest
export integerCheck
export createBranch

end
