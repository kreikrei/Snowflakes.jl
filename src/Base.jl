using XLSX
using DataFrames

function base(path::String)
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

    K = Dict{Int64,veh}() #INITIATE VEHICLES

    for k in eachrow(data[:vehicles]) #ITERATE OVER DATA
        K[k.id] = veh(
            k.name , k.type ,
            parse.(Int64,split(k.cover)) , parse.(Int64,split(k.loadp)) , k.freq , k.Q ,
            k.varq , k.vardq , k.vard , k.fix
        )
    end

    dt = (V=V,K=K) #WRAPPING RESULT

    return dt
end

function report(res::NamedTuple)
    uniqueVtx = unique([res.V[i].type for i in keys(res.V)])
    uniqueVeh = unique([res.K[k].type for k in keys(res.K)])

    vtxType = DataFrame()

    stats = status(
        length(res.V) ,
        length(res.K) ,
        uniqueVtx ,
        uniqueVeh
        #vehicle_type_breakdown
        #vertice_type_breakdown
        #unique_cover
        #unique_loadp
    )

    return stats
end
