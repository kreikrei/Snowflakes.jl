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
    y = r.y.data .- p.y.data
    z = r.z.data .- p.z.data
    x = r.x.data .- p.x.data

    yB = isempty(filter(p -> p < 0, y))
    zB = isempty(filter(p -> p < 0, z))
    xB = isempty(filter(p -> p < 0, x))

    if yB && zB && xB
        return true
    else
        return false
    end
end

function Q(vector::col,idx::NamedTuple;R::Dict{Int64,col})
    set = Vector{Int64}()

    for r in keys(R)
        if dominance(qvec(r,idx.k,idx.t;R=R),vector)
            push!(set,r)
        end
    end

    return set
end

function positiveComp(q::col)
    return sum(q.y) + sum(q.z) + sum(q.x)
end

function separate(n::node)
    collection = DataFrame(
        q = col[], idx = NamedTuple[], val = Float64[],
        set = Vector{Int64}[],
        positive=Float64[]
    )

    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    for r in keys(R), k in keys(b().K), t in b().T
        tot = sum(θ[q,k,t] for q in Q(qvec(r,k,t;R=R),(k=k,t=t);R=R))
        if !(abs(round(tot) - tot) < 1e-15) #isinteger(tot)
            append!(collection,
                DataFrame(
                    q = qvec(r,k,t;R=R),
                    idx = (k=k,t=t),
                    val = tot,
                    set = [Q(qvec(r,k,t;R=R),(k=k,t=t);R=R)],
                    positive = positiveComp(qvec(r,k,t;R=R))
                )
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
                            br, s.idx, s.q,
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
