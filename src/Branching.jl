# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function qvec(r::Int64,k::Int64,t::Int64;R::Dict{Int64,col})
    return col(
        R[r].u[:,k,t],
        R[r].v[:,k,t],
        R[r].l[:,:,k,t],
        R[r].y[:,k,t],
        R[r].z[:,k,t],
        R[r].x[:,:,k,t]
    )
end

function dominance(r::col,p::col)
    u = r.u.data .- p.u.data
    v = r.v.data .- p.v.data
    l = r.l.data .- p.l.data
    y = r.y.data .- p.y.data
    z = r.z.data .- p.z.data
    x = r.x.data .- p.x.data

    uB = isempty(filter(p -> p < 0, u))
    vB = isempty(filter(p -> p < 0, v))
    lB = isempty(filter(p -> p < 0, l))
    yB = isempty(filter(p -> p < 0, y))
    zB = isempty(filter(p -> p < 0, z))
    xB = isempty(filter(p -> p < 0, x))

    if uB && vB && lB && yB && zB && xB
        return true
    else
        return false
    end
end

function Q(q::col;R::Dict{Int64,col})
    set = Vector{NamedTuple}()

    for r in keys(R), k in keys(b().K), t in b().T
        if dominance(qvec(r,k,t;R=R),q)
            push!(set,(r=r,k=k,t=t))
        end
    end

    return set
end

function positiveComp(q::col)
    return sum(q.u) + sum(q.v) + sum(q.l) + sum(q.y) + sum(q.z) + sum(q.x)
end

function separate(n::node)
    collection = DataFrame(
        q = col[], val = Float64[],
        set = Vector{NamedTuple}[],
        positive=Float64[]
    )

    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    for r in keys(R), k in keys(b().K), t in b().T
        tot = sum(θ[q.r,q.k,q.t] for q in Q(qvec(r,k,t;R=R);R=R))
        if !(abs(round(tot) - tot) < 1e-15) #isinteger(tot)
            append!(collection,
                DataFrame(
                    q = qvec(r,k,t;R=R),
                    val = tot,
                    set = [Q(qvec(r,k,t;R=R);R=R)],
                    positive = positiveComp(qvec(r,k,t;R=R)))
            )
        end
    end

    return collection
end

function integerCheck(n::node)
    integer = true
    sol = origin(n)

    for i in keys(b().V), k in keys(b().K), t in b().T
        if !isinteger(sol.z[i,k,t]) || !isinteger(sol.y[i,k,t])
            integer = false
            break
        end
    end

    return integer
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
                            br, s.q,
                            if br == "upper"
                                floor(s.val)
                            else
                                ceil(s.val)
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
