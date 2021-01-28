# =========================================================================
#    BASIC STRUCTURES
# =========================================================================

struct vtx
    #IDENTIFIERS
    name::String
    type::String

    #VALS
    x::Float64 #xcoor
    y::Float64 #ycoor
    MAX::Float64 #max inventory level
    MIN::Float64 #min inventory level
    START::Float64 #starting inventory level

    #COSTS
    h::Float64 #inv cost per unit
end

struct veh
    #IDENTIFIERS
    name::String
    type::String

    #CHAR
    cover::Vector{Int64}
    freq::Int64
    Q::Int64

    #COSTS
    vx::Float64
    vl::Float64
    fp::Float64
    fd::Float64
end

struct dt
    V::Dict{Int64,vtx}
    dist::JuMP.Containers.DenseAxisArray
    K::Dict{Int64,veh}
    T::Vector{Int64}
    d::JuMP.Containers.DenseAxisArray
end

struct col
    #q related
    q::JuMP.Containers.SparseAxisArray
    u::JuMP.Containers.SparseAxisArray
    v::JuMP.Containers.SparseAxisArray
    l::JuMP.Containers.SparseAxisArray

    #0-1 related
    p::JuMP.Containers.SparseAxisArray
    y::JuMP.Containers.SparseAxisArray
    z::JuMP.Containers.SparseAxisArray
    x::JuMP.Containers.SparseAxisArray
end

struct bound
    idx::NamedTuple
    val::Int64
end

struct dval
    λ::JuMP.Containers.DenseAxisArray
    γ::JuMP.Containers.DenseAxisArray
    δ::JuMP.Containers.SparseAxisArray
end

struct stabilizer
    slCoeff::Float64
    suCoeff::Float64
    slLim::JuMP.Containers.DenseAxisArray
    suLim::JuMP.Containers.DenseAxisArray
end

struct node
    #IDENTIFIER
    parent::UUID
    self::UUID

    #Dynamaic SET
    bounds::Vector{bound}
    columns::Vector{col}

    #SUPPORT
    stab::stabilizer

    #STATUS
    status::Vector{String}
end
