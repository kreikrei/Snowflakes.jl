# =========================================================================
#    BRANCHING MECHANISMS
# =========================================================================

function separate(n::Snowflakes.node)
    #PROTOTYPE SEPARATION
    R = Dict(1:length(n.columns) .=> n.columns)
    θ = value.(master(n).obj_dict[:θ])

    F = Vector{NamedTuple}()
    for r in keys(R), k in keys(b().K), t in b().T
        if θ[r,k,t] - floor(θ[r,k,t]) > 0
            push!(F,(r=r,k=k,t=t))
        end
    end #semua q yang nilai θ-nya masih fractional

    qF = DataFrame(q = Symbol[], i = Int64[], val = Int64[])
    for f in F, q in [:y,:z,:v,:u,], i in keys(b().V)
        append!(qF,
            DataFrame(
                q = q,
                i = i,
                val = getproperty(R[f.r],q)[i,f.k,f.t]
            )
        )
    end #collect existing values of the variable (q,i)

    return qF
end

function vtest(q::Symbol,i::Int64,qF::DataFrame)
    test_v = Vector{Int64}()

    distinct = filter(p -> p.i == i && p.q == q,qF).val
    sort!(distinct)

    for i in 1:length(distinct)-1
        comp = ceil((distinct[i] + distinct[i+1]) / 2)
        push!(test_v, comp)
    end

    res = Vector{β}()

    for v in unique!(test_v)
        push!(res,β(q,i,"≳",v))
    end

    return res
end

function integerCheck(n::node)
    integer = true
    sol = origin(n)

    for i in keys(b().V), k in keys(b().K), t in b().T
        if !isinteger(sol.y[i,k,t])
            integer = false
            break
        end
    end

    return integer
end

function createBranch(n::node)
    branches = Vector{node}()
    seeds = separate(n)


    for br in ["≳","≲"]
        push!(branches,
            node(
                n.self, #parent
                uuid1(), #self

                vcat(n.bounds, #bounds
                    bound(
                        br, first(seeds),
                        if br == "≲"
                            floor(last(seeds))
                        else
                            ceil(last(seeds))
                        end
                    )
                ),
                n.columns, #columns
                initStab(), #stabilizer
                ["UNVISITED"]
            )
        )
    end

    return branches
end
