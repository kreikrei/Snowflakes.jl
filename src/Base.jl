# =========================================================================
#    BASIC DATA GENERATION
# =========================================================================

const template_data = Ref{Any}(nothing)
extract() = template_data[] #buat manggil dt default

function base!(path::String) #extract from excel
    xf = XLSX.readxlsx(path) #READ WORKSHEET
    data = Dict{Symbol,DataFrame}() #DATAFRAME DICT

    for sheets in XLSX.sheetnames(xf) #TURN SHEETS INTO DATAFRAME
        df = DataFrame(XLSX.gettable(xf[sheets])...) #TRANSFORM TABLE INTO DATAFRAME
        data[Symbol(sheets)] = df #DEFINE THE NAME FROM THE WORKSHEET
    end

    V = Dict{Int64,vtx}() #INITIATE VERTICES
    for v in eachrow(data[:vertices]) #ITERATE OVER DATA
        V[v.id] = vtx(
            v.name, v.type,
            v.x, v.y, v.MAX, v.MIN, v.START,
            v.h
        )
    end

    dist = JuMP.Containers.DenseAxisArray{Float64}(undef, keys(V), keys(V))
    for i in keys(V), j in keys(V)
        if i != j
            dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],6378.137)
        else
            dist[i,j] = 999999999
        end
    end

    K = Dict{Int64,veh}() #INITIATE VEHICLES
    for k in eachrow(data[:vehicles]) #ITERATE OVER DATA
        K[k.id] = veh(
            k.name, k.type,
            parse.(Int64,split(k.cover)), parse.(Int64,split(k.loadp)), k.freq, k.Q,
            k.vx, k.vl, k.fp, k.fd
        )
    end

    T = collect( #range from starting month for duration
        range(
            last(data[:periods].start),
            length = last(data[:periods].T),
            step = 1
        )
    )

    d = JuMP.Containers.DenseAxisArray(
        Array{Float64}(data[:demands][:,string.(T)]), #dataset
        Array{Int64}(data[:demands].point), #dims 1
        T #dims 2
    )

    res = dt(V,dist,K,T,d) #WRAPPING RESULT TO TYPE dt

    return template_data[] = res
end

function stats(res = extract())
    uniqueVtx = unique([res.V[i].type for i in keys(res.V)])
    uniqueVeh = unique([res.K[k].type for k in keys(res.K)])

    vtxType = Dict{String,Vector{Int64}}()
    for t in uniqueVtx
        vtxType[t] = sort(collect(keys(filter(p -> last(p).type == t , res.V))))
    end
    vehType = Dict{String,Vector{Int64}}()
    for t in uniqueVeh
        vehType[t] = sort(collect(keys(filter(p -> last(p).type == t , res.K))))
    end

    loadp = Vector{Int64}()
    cover = Vector{Int64}()
    for k in keys(res.K)
        for l in res.K[k].loadp
            push!(loadp,l)
        end
        for c in res.K[k].cover
            push!(cover,c)
        end

        unique!(loadp)
        unique!(cover)
    end

    stats = (
        number_of_vertices = length(res.V),
        number_of_vehicles = length(res.K),
        unique_types_vtx = uniqueVtx,
        unique_types_veh = uniqueVeh,
        type_breakdown_vtx = vtxType,
        type_breakdown_veh = vehType,
        cover_list = cover,
        loadp_list = loadp
    )

    return stats
end

function initStab(res = extract();slC::Float64,suC::Float64)
    slackCoeff = slC
    surpCoeff = suC
    slackLim = abs.(res.d)
    surpLim = abs.(res.d)

    return stabilizer(slackCoeff,surpCoeff,slackLim,surpLim)
end

function root(res = extract();slC::Float64,suC::Float64)
    id = uuid1()
    root = node(
        id, id,
        res,
        Vector{bound}(),Vector{col}(),
        initStab(res,slC = slC, suC = suC),
        ["UNVISITED"]
    )

    return root
end
