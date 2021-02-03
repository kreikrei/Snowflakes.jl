# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function integerCheck(n::node)
    integer = true

    sol = origin(n)

    #z check
    for z in sol.z, y in sol.y
        if !isinteger(z)
            integer = false
            break
        end

        if !isinteger(y)
            integer = false
            break
        end
    end

    return integer
end

function fractionalPart(var::Float64)
    a = var - floor(var)
    b = ceil(var) - var

    return min(a,b)
end

function exploreFrac(n::node)
    R = Dict(1:length(n.columns) .=> n.columns)
    mp = master(n)
    θ = value.(mp.obj_dict[:θ])
    collection = DataFrame(
        var=Symbol[],
        idx=NamedTuple[],
        e=Int64[],
        r=Vector{Int64}[],
        val=Float64[],
        frac=Float64[]
    )

    for i in keys(b().V), k in keys(b().K), t in b().T
        for e in collect(1:b().K[k].freq)
            iter = keys(filter(m -> last(m).y[i,k,t] >= e,R)) #iter y

            if !isempty(iter)
                tot = sum(
                    θ[r,k,t]
                    for r in iter
                )

                if tot > 0 && !isinteger(tot)
                    append!(collection,DataFrame(
                            var=:y, idx=(i=i,k=k,t=t),
                            e=e, r=[collect(iter)], val=tot,
                            frac=fractionalPart(tot)
                        )
                    )
                end
            end

            iter = keys(filter(m -> last(m).z[i,k,t] >= e,R)) #iter z

            if !isempty(iter)
                tot = sum(
                    θ[r,k,t]
                    for r in iter
                )

                if tot > 0 && !isinteger(tot)
                    append!(collection,DataFrame(
                            var=:z, idx=(i=i,k=k,t=t),
                            e=e, r=[collect(iter)], val=tot,
                            frac=fractionalPart(tot)
                        )
                    )
                end
            end
        end
    end

    return collection
end

function createBranch(n::node,seeds::DataFrame)
    branches = Vector{node}()

    for s in eachrow(seeds)
        for br in ["upper","lower"]
            push!(branches,
                node(
                    n.self, #parent
                    uuid1(), #self
                    vcat(n.bounds, #bounds
                        bound(
                            s.var, br,
                            s.idx, s.e,
                            if br == "upper"
                                ceil(s.val)
                            else
                                floor(s.val)
                            end
                        )
                    ),
                    n.columns, #columns
                    initStab(), #stabilizer
                    ["UNVISITED"]
                )
            )
        end
    end

    return branches
end
