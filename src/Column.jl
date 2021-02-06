# =========================================================================
#    COLUMN GENERATION MECHANISMS
# =========================================================================

function Q(key,R::Dict)
    if isa(key,Tuple)
        q = keys(filter(p -> first(last(p)) == key,R))
    elseif isa(key,β)
        q = keys(filter(p -> dominance(last(last(p)),key.idx,key.value),R))
    elseif isa(key,Vector{β})
        q = collect(keys(R))
        for b in key
            intersect!(q,Q(b,R))
        end
    end

    return q
end

function dominance(p::col,i::Int64,val::Int64)
    if p.y[i] >= val
        return true
    else
        return false
    end
end

function f(B,R,θ)
    if isempty(B)
        return sum(θ[r,first(R[r].first),last(R[r].first)] -
        floor(θ[r,first(R[r].first),last(R[r].first)])
            for r in keys(R)
        )
    else
        if !isempty(Q(B,R))
            return sum(θ[r,first(R[r].first),last(R[r].first)] -
                floor(θ[r,first(R[r].first),last(R[r].first)])
                for r in Q(B,R))
        else
            return 0
        end
    end
end

function master(n::node)
    mp = Model(get_default_optimizer())
    if silent()
        set_silent(mp)
    end

    R = Dict(1:length(n.columns) .=> n.columns)
    B = Dict(1:length(n.bounds) .=> n.bounds)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    @variable(mp, θ[keys(R), keys(b().K), b().T] >= 0)
    @variable(mp, I[keys(b().V), vcat(first(b().T) - 1, b().T)])
    @variable(mp, 0 <= slack[i = keys(b().V), t = b().T] <= n.stab.slLim[i,t])
    @variable(mp, 0 <= surp[i = keys(b().V), t = b().T] <= n.stab.suLim[i,t])

    @objective(mp, Min,
        sum(
            θ[r,k,t] * sum(
                sum(
                    b().dist[i,j] * (
                        b().K[k].vx * last(R[r]).x[i,j] +
                        b().K[k].vl * last(R[r]).l[i,j]
                    )
                    for i in b().K[k].cover, j in b().K[k].cover
                ) +
                sum(
                    b().K[k].fd * last(R[r]).u[i]
                    for i in b().K[k].cover
                ) +
                sum(
                    b().K[k].fp * last(R[r]).z[i]
                    for i in b().K[k].cover
                )
                for r in Q((k,t),R)
            )
            for k in keys(b().K), t in b().T
        ) + #column costs
        sum(
            b().V[i].h * I[i,t]
            for i in keys(b().V), t in b().T
        ) + #inventory costs
        sum(
            n.stab.slCoeff * slack[i,t]
            for i in keys(b().V), t in b().T
        ) - #stabilizer
        sum(
            n.stab.suCoeff * surp[i,t]
            for i in keys(b().V), t in b().T
        ) #stabilizer
    )

    @constraint(mp, λ[i = keys(b().V), t = b().T],
        I[i,t - 1] + sum(last(R[r]).u[i] * θ[r,k,t] for k in keys(b().K), r in Q((k,t),R))
        + slack[i,t] - surp[i,t] ==
        sum(last(R[r]).v[i] * θ[r,k,t] for k in keys(b().K), r in Q((k,t),R)) +
        b().d[i,t] + I[i,t]
    )

    @constraint(mp, δ[k = keys(b().K), t = b().T],
        sum(θ[r,k,t] for r in Q((k,t),R)) <= 1 #convexity constraint
    )

    @constraint(mp, [i = keys(b().V), t = b().T],
        b().V[i].MIN <= I[i,t] <= b().V[i].MAX #inventory capacity interval
    )

    @constraint(mp, [i = keys(b().V)],
        I[i,first(b().T)-1] == b().V[i].START #starting inventory level
    )

    # ================================
    #    BOUND GENERATOR
    # ================================
    ≲ = filter(b -> last(b).type == "≲",B)
    ≳ = filter(b -> last(b).type == "≳",B)

    @constraint(mp, μ[b = keys(≲)],
        sum(θ[q,first(R[r].first),last(R[r].first)] for q in Q(B[b].B,R)) <= B[b].κ
    )

    @constraint(mp, ν[b = keys(≳)],
        sum(θ[q,first(R[r].first),last(R[r].first)] for q in Q(B[b].B,R)) >= B[b].κ
    )

    optimize!(mp)

    return mp
end

function getDuals(mp::Model)
    λ = dual.(mp.obj_dict[:λ])
    δ = dual.(mp.obj_dict[:δ])
    μ = dual.(mp.obj_dict[:μ])
    ν = dual.(mp.obj_dict[:ν])

    return dval(λ,δ,μ,ν)
end

function sub(n::node,duals::dval)
    sp = Model(get_default_optimizer())
    if silent()
        set_silent(sp)
    end

    R = Dict(1:length(n.columns) .=> n.columns)
    B = Dict(1:length(n.bounds) .=> n.bounds)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    q = Dict{Tuple,Snowflakes.col}()

    for k in keys(b().K), t in b().T
        q[(k,t)] = Snowflakes.col(
            @variable(sp, [i = keys(b().V)], Int), #u
            @variable(sp, [i = keys(b().V)], Int), #v
            @variable(sp, [i = collect(keys(b().V)), j = collect(keys(b().V))], Int), #l
            @variable(sp, [i = keys(b().V)], Int), #y
            @variable(sp, [i = keys(b().V)], Int), #z
            @variable(sp, [i = collect(keys(b().V)), j = collect(keys(b().V))], Int) #x
        )

        @constraint(sp, [i = keys(b().V)], 0 <= q[(k,t)].u[i] <=
            b().K[k].freq * length(b().K[k].cover) * b().K[k].Q
        )
        @constraint(sp, [i = keys(b().V)], 0 <= q[(k,t)].v[i] <=
            b().K[k].freq * b().K[k].Q
        )
        @constraint(sp, [i = collect(keys(b().V)), j = collect(keys(b().V))], 0 <=
            q[(k,t)].l[i,j] <= b().K[k].freq * length(b().K[k].cover) * b().K[k].Q
        )
        @constraint(sp, [i = keys(b().V)], 0 <= q[(k,t)].y[i] <=
            b().K[k].freq * length(b().K[k].cover)
        )
        @constraint(sp, [i = keys(b().V)], 0 <= q[(k,t)].z[i] <=
            b().K[k].freq
        )
        @constraint(sp, [i = collect(keys(b().V)), j = collect(keys(b().V))], 0 <=
            q[(k,t)].x[i,j] <= b().K[k].freq * length(b().K[k].cover)
        )
    end

    @variable(sp, o[keys(B)], Bin)

    @objective(sp, Min,
        sum(
            sum(
                b().dist[i,j] * (
                    b().K[k].vx * q[(k,t)].x[i,j] +
                    b().K[k].vl * q[(k,t)].l[i,j]
                )
                for i in b().K[k].cover, j in b().K[k].cover
            ) +
            sum(
                b().K[k].fd * q[(k,t)].u[i]
                for i in b().K[k].cover
            ) +
            sum(
                b().K[k].fp * q[(k,t)].z[i]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                (q[(k,t)].u[i] - q[(k,t)].v[i]) * duals.λ[i,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            duals.δ[k,t]
            for k in keys(b().K), t in b().T
        ) -
        sum(
            o[b] * duals.μ[b]
            for b in keys(filter(b -> last(b).type == "upper",B))
        ) -
        sum(
            o[b] * duals.ν[b]
            for b in keys(filter(b -> last(b).type == "lower",B))
        )
    )

    @constraint(sp, [k = keys(b().K), t = b().T],
        sum(q[(k,t)].u[i] - q[(k,t)].v[i] for i in b().K[k].cover) == 0 #all pickup delivered
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        q[(k,t)].u[i] <= b().K[k].Q * q[(k,t)].y[i] #u-y corr
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        q[(k,t)].v[i] <= b().K[k].Q * q[(k,t)].z[i] #v-z corr
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(q[(k,t)].x[j,i] for j in b().K[k].cover) == q[(k,t)].y[i] + q[(k,t)].z[i] #traverse in
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(q[(k,t)].x[i,j] for j in b().K[k].cover) == q[(k,t)].y[i] + q[(k,t)].z[i] #traverse out
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(q[(k,t)].l[j,i] for j in b().K[k].cover) -
        sum(q[(k,t)].l[i,j] for j in b().K[k].cover) == q[(k,t)].u[i] - q[(k,t)].v[i] #load balance
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, j = b().K[k].cover, t = b().T],
        q[(k,t)].l[i,j] <= b().K[k].Q * q[(k,t)].x[i,j] #l-x corr
    )

    @constraint(sp, [k = keys(b().K),
        n = [i for i in keys(b().V) if !(i in b().K[k].cover)], t = b().T],
        q[(k,t)].z[n] == 0
    )

    @constraint(sp, [k = keys(b().K),
        n = [i for i in keys(b().V) if !(i in b().K[k].cover)], t = b().T],
        q[(k,t)].y[n] == 0
    )

    # ================================
    #    BOUND GENERATOR
    # ================================
    for id in keys(B)
        y_sup = filter(p -> B[id].vector.y[p] > 0,collect(B[id].vector.y.axes[1])) #set of i
        𝕪 = @variable(sp, [y_sup], Bin)
        @constraint(sp, [i = y_sup],
            B[id].vector.y[i] * 𝕪[i] <= y[i,B[id].idx.k,B[id].idx.t]
        )
        @constraint(sp, [i = y_sup],
            y[i,B[id].idx.k,B[id].idx.t] <= (B[id].vector.y[i] - 1) +
            (length(b().K[B[id].idx.k].cover) * b().K[B[id].idx.k].freq -
            B[id].vector.y[i] + 1) * 𝕪[i]
        )
        @constraint(sp, [i = y_sup], o[id] <= 𝕪[i])

        z_sup = filter(p -> B[id].vector.z[p] > 0,collect(B[id].vector.z.axes[1]))
        𝕫 = @variable(sp, [z_sup], Bin)
        @constraint(sp, [i = z_sup],
            B[id].vector.z[i] * 𝕫[i] <= z[i,B[id].idx.k,B[id].idx.t]
        )
        @constraint(sp, [i = z_sup],
            z[i,B[id].idx.k,B[id].idx.t] <= (B[id].vector.z[i] - 1) +
            (b().K[B[id].idx.k].freq -
            B[id].vector.z[i] + 1) * 𝕫[i]
        )
        @constraint(sp, [i = z_sup], o[id] <= 𝕫[i])

        x_sup = filter(p -> B[id].vector.x[first(p),last(p)] > 0,
        collect(combinations(B[id].vector.x.axes[1],2)))
        𝕩 = @variable(sp, [x_sup], Bin)
        @constraint(sp, [i = x_sup],
            B[id].vector.x[first(i),last(i)] * 𝕩[i] <=
            x[first(i),last(i),B[id].idx.k,B[id].idx.t]
        )
        @constraint(sp, [i = x_sup],
            x[first(i),last(i),B[id].idx.k,B[id].idx.t] <=
            (B[id].vector.x[first(i),last(i)] - 1) +
            (length(b().K[B[id].idx.k].cover) * b().K[B[id].idx.k].freq -
            B[id].vector.x[first(i),last(i)] + 1) * 𝕩[i]
        )
        @constraint(sp, [i = x_sup], o[id] <= 𝕩[i])

        @constraint(sp, o[id] >= 1 - (
                sum(1 - 𝕪[i] for i in y_sup) +
                sum(1 - 𝕫[i] for i in z_sup) +
                sum(1 - 𝕩[i] for i in x_sup)
            )
        )
    end

    optimize!(sp)

    return sp,q
end

function getCols(sp::Dict{Tuple,Snowflakes.col})
    q = Dict{Tuple,Snowflakes.col}()

    for k in keys(b().K), t in b().T
        q[(k,t)] = Snowflakes.col(
            value.(sp[(k,t)].u),
            value.(sp[(k,t)].v),
            value.(sp[(k,t)].l),
            value.(sp[(k,t)].y),
            value.(sp[(k,t)].z),
            value.(sp[(k,t)].x)
        )
    end

    return [i for i in q]
end

function updateStab!(stab::stabilizer,param::Float64)
    for i in first(stab.slLim.axes),t in last(stab.slLim.axes)
        stab.slLim[i,t] = param * stab.slLim[i,t]
        if stab.slLim[i,t] < 1
            stab.slLim[i,t] = 0
        end
    end

    for i in first(stab.suLim.axes),t in last(stab.suLim.axes)
        stab.suLim[i,t] = param * stab.suLim[i,t]
        if stab.suLim[i,t] < 1
            stab.suLim[i,t] = 0
        end
    end

    return stab
end

function checkStab(mp::Model)
    s = sum(value.(mp.obj_dict[:slack])) + sum(value.(mp.obj_dict[:surp]))

    return s
end

function colGen(n::node;maxCG::Float64,track::Bool)
    terminate = false
    iter = 0
    mem = 0

    while !terminate
        if iter < maxCG
            mp = master(n)

            if has_values(mp) && has_duals(mp)
                if track #print master problem obj
                    println("obj: $(objective_value(mp))")
                end

                duals = getDuals(mp)
                sp = sub(n,duals)

                if track #print subproblem price
                    println("price: $(objective_value(sp[1]))")
                end

                if (isapprox(objective_value(sp[1]),0,atol = 1e-8) ||
                    objective_value(sp[1]) > 0)
                    if isapprox(checkStab(mp),0,atol = 1e-8)
                        terminate = true #action
                        push!(n.status,"EVALUATED") #report
                        if track
                            println("EVALUATED")
                        end
                    else
                        updateStab!(n.stab,0.5) #action
                        push!(n.status,"STABILIZED") #report
                        if track
                            println("STABILIZED")
                        end
                    end
                else
                    append!(n.columns,getCols(sp[2])) #action
                    push!(n.status,"ADD_COLUMN") #report
                    if track
                        println("ADD_COLUMN")
                    end
                end

                iter += 1 #iteration update
            else
                terminate = true #action
                push!(n.status,"NO_SOLUTION")
                if track
                    println("NO_SOLUTION")
                end
            end
        else
            terminate = true #action
            push!(n.status,"EVALUATED") #report
            if track
                println("EVALUATED")
            end
        end
    end

    if n.status[end] == "NO_SOLUTION"
        println("NODE $(n.self) FAILED.")
    else
        println("NODE $(n.self) FINISHED.")
    end

    return n
end
