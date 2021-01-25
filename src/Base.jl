# ==============================================================================
#    BASIC DATA GENERATION
# ==============================================================================

function base(path::String) #extract from excel
    xf = XLSX.readxlsx(path) #READ WORKSHEET
    data = Dict{Symbol,DataFrame}() #DATAFRAME DICT

    for sheets in XLSX.sheetnames(xf) #TURN SHEETS INTO DATAFRAME
        df = DataFrame(XLSX.gettable(xf[sheets])...) #TRANSFORM TABLE INTO DATAFRAME
        data[Symbol(sheets)] = df #DEFINE THE NAME FROM THE WORKSHEET
    end

    V = Dict{Int64,vtx}() #INITIATE VERTICES
    for v in eachrow(data[:vertices]) #ITERATE OVER DATA
        V[v.id] = vtx(
            v.name , v.type,
            v.x , v.y , v.MAX , v.MIN , v.START ,
            v.h
        )
    end

    dist = distance(V)

    K = Dict{Int64,veh}() #INITIATE VEHICLES
    for k in eachrow(data[:vehicles]) #ITERATE OVER DATA
        K[k.id] = veh(
            k.name , k.type ,
            parse.(Int64,split(k.cover)) , parse.(Int64,split(k.loadp)) , k.freq , k.Q ,
            k.varq , k.vardq , k.vard , k.fix
        )
    end

    T = collect(
        range(
            last(data[:periods].start) ,
            length = last(data[:periods].T) ,
            step = 1
        )
    ) #range starting from starting month for duration

    d = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),T)
    for i in keys(V)
        row = @from x in data[:demands] begin
            @where x.point == i
            @select x
            @collect DataFrame
        end

        for t in T
            d[i,t] = row[1,2:ncol(row)-3][t]
        end
    end

    dt = (V=V,K=K,T=T,d=d,dist=dist) #WRAPPING RESULT

    return dt
end

function report(res::NamedTuple)
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

    stats = status(
        length(res.V) ,
        length(res.K) ,
        uniqueVtx ,
        uniqueVeh ,
        vtxType ,
        vehType ,
        cover ,
        loadp
    )

    return stats
end

function distance(V::Dict)
    dist = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),collect(keys(V)))
    for i in keys(V), j in keys(V)
        if i != j
            dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],6378.137)
        else
            dist[i,j] = 999999999
        end
    end

    return dist
end
