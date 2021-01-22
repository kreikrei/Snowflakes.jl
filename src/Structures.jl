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
    loadp::Vector{Int64}
    freq::Int64
    Q::Int64

    #COSTS
    varq::Float64
    vardq::Float64
    vard::Float64
    fix::Float64
end

struct status
    #FOR ASSESSING A DATASET
    number_of_vertices::Int64
    number_of_vehicles::Int64
    unique_vertices_type::Vector{String}
    unique_vehicle_types::Vector{String}
    vertex_type_breakdown::Dict{String,Vector{Int64}}
    vehicle_type_breakdown::Dict{String,Vector{Int64}}
    unique_cover::Vector{Int64}
    unique_loadp::Vector{Int64}
end
