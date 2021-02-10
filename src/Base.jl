# =========================================================================
#    BASIC DATA GENERATION
# =========================================================================

const template_data = Ref{Any}(nothing)
b() = template_data[] #buat manggil dt default

const max_i = Ref{Any}(nothing)
qmax() = max_i[] #buat manggil nilai max di i

function imax(res=b())
    val = col(
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(res.V)), #u
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(res.V)), #v
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(res.V), keys(res.V)), #empty l
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(res.V)), #y
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(res.V)), #z
        JuMP.Containers.DenseAxisArray{Float64}(undef, keys(res.V), keys(res.V)) #empty x
    )

    for i in keys(res.V)
        u = [length(res.K[k].cover) * res.K[k].freq * res.K[k].Q for k in keys(res.K)
        if i in res.K[k].cover]
        v = [res.K[k].freq * res.K[k].Q for k in keys(res.K) if i in res.K[k].cover]

        y = [length(res.K[k].cover) * res.K[k].freq for k in keys(res.K)
        if i in res.K[k].cover]
        z = [res.K[k].freq for k in keys(res.K) if i in res.K[k].cover]

        val.u[i] = findmax(u)[1]
        val.v[i] = findmax(v)[1]
        
        val.y[i] = findmax(y)[1]
        val.z[i] = findmax(z)[1]
    end

    return max_i[] = val
end

function extract!(path::String) #extract from excel
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
            parse.(Int64,split(k.cover)), k.freq, k.Q,
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

function stats(res = b())
    uniqueVtx = unique([res.V[i].type for i in keys(res.V)])
    uniqueVeh = unique([res.K[k].type for k in keys(res.K)])

    vtxType = Dict{String,Vector{Int64}}()
    vehType = Dict{String,Vector{Int64}}()
    for t in uniqueVtx
        vtxType[t] = sort(collect(keys(filter(p -> last(p).type == t , res.V))))
    end
    for t in uniqueVeh
        vehType[t] = sort(collect(keys(filter(p -> last(p).type == t , res.K))))
    end

    cover = Vector{Int64}()
    for k in keys(res.K)
        for c in res.K[k].cover
            push!(cover,c)
        end
        unique!(cover)
    end

    dems = [mean(res.d[:,t]) for t in res.T]

    return status(
        length(res.V),
        length(res.K),
        uniqueVtx,
        uniqueVeh,
        vtxType,
        vehType,
        cover,
        dems
    )
end

function initStab(res = b())
    slackCoeff = sl_C()
    surpCoeff = su_C()
    slackLim = abs.(res.d)
    surpLim = abs.(res.d)

    return stabilizer(slackCoeff,surpCoeff,slackLim,surpLim)
end

function root(res = b())
    id = uuid1()
    root = node(
        id, id,
        Vector{bound}(),Vector{col}(),
        initStab(),
        ["UNVISITED"]
    )

    return root
end

function origin(n::node)
    z = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    y = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    u = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    v = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(b().V),keys(b().K),b().T)
    x = JuMP.Containers.DenseAxisArray{Float64}(
        undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K)),b().T
    )
    l = JuMP.Containers.DenseAxisArray{Float64}(
        undef,collect(keys(b().V)),collect(keys(b().V)),collect(keys(b().K)),b().T
    )

    R = Dict(1:length(n.columns) .=> n.columns)
    mp = master(n)
    θ = mp.obj_dict[:θ]

    for i in keys(b().V), k in keys(b().K), t in b().T
        z[i,k,t] = value(sum(R[r].z[i,k,t] * θ[r,k,t] for r in keys(R)))
        y[i,k,t] = value(sum(R[r].y[i,k,t] * θ[r,k,t] for r in keys(R)))
        u[i,k,t] = value(sum(R[r].u[i,k,t] * θ[r,k,t] for r in keys(R)))
        v[i,k,t] = value(sum(R[r].v[i,k,t] * θ[r,k,t] for r in keys(R)))

        for j in keys(b().V)
            x[i,j,k,t] = value(sum(R[r].x[i,j,k,t] * θ[r,k,t] for r in keys(R)))
            l[i,j,k,t] = value(sum(R[r].l[i,j,k,t] * θ[r,k,t] for r in keys(R)))
        end
    end

    return col(
        u,v,l,
        y,z,x
    )
end
