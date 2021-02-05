# =========================================================================
#    COLUMN GENERATION MECHANISMS
# =========================================================================

function master(n::node)
    mp = buildMaster(n)
    optimize!(mp)

    return mp
end

function buildMaster(n::node)
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
            θ[r,k,t] * (
                sum(
                    b().dist[i,j] * (
                        b().K[k].vx * R[r].x[i,j,k,t] +
                        b().K[k].vl * R[r].l[i,j,k,t]
                    )
                    for i in b().K[k].cover, j in b().K[k].cover
                ) +
                sum(
                    b().K[k].fd * R[r].u[i,k,t]
                    for i in b().K[k].cover
                ) +
                sum(
                    b().K[k].fp * R[r].z[i,k,t]
                    for i in b().K[k].cover
                )
            )
            for r in keys(R), k in keys(b().K), t in b().T
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
        I[i,t - 1] + sum(R[r].u[i,k,t] * θ[r,k,t] for r in keys(R), k in keys(b().K))
        + slack[i,t] - surp[i,t] ==
        sum(R[r].v[i,k,t] * θ[r,k,t] for r in keys(R), k in keys(b().K)) +
        b().d[i,t] + I[i,t]
    )

    @constraint(mp, δ[k = keys(b().K), t = b().T],
        sum(θ[r,k,t] for r in keys(R)) <= 1 #convexity constraint
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
    uB = filter(b -> last(b).type == "upper",B)
    lB = filter(b -> last(b).type == "lower",B)

    @constraint(mp, μ[b = keys(uB)],
        sum(θ[q,B[b].idx.k,B[b].idx.t] for q in Q(B[b].vector,B[b].idx;R=R)) <= B[b].value
    )

    @constraint(mp, ν[b = keys(lB)],
        sum(θ[q,B[b].idx.k,B[b].idx.t] for q in Q(B[b].vector,B[b].idx;R=R)) >= B[b].value
    )

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
    sp = buildSub(n,duals)
    optimize!(sp)

    return sp
end

function buildSub(n::node,duals::dval)
    sp = Model(get_default_optimizer())
    if silent()
        set_silent(sp)
    end

    R = Dict(1:length(n.columns) .=> n.columns)
    B = Dict(1:length(n.bounds) .=> n.bounds)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    @variable(sp, u[keys(b().V), keys(b().K), b().T] >= 0, Int)
    @variable(sp, v[keys(b().V), keys(b().K), b().T] >= 0, Int)
    @variable(sp,
        l[collect(keys(b().V)), collect(keys(b().V)), collect(keys(b().K)), b().T] >= 0, Int
    )

    @variable(sp, y[keys(b().V), keys(b().K), b().T] >= 0, Int)
    @variable(sp, 0 <= z[i = keys(b().V), k = keys(b().K), t = b().T] <= b().K[k].freq, Int)
    @variable(sp,
        x[collect(keys(b().V)), collect(keys(b().V)), collect(keys(b().K)), b().T] >= 0, Int
    )

    @variable(sp, o[keys(B)], Bin)

    @objective(sp, Min,
        sum(
            sum(
                b().dist[i,j] * (
                    b().K[k].vx * x[i,j,k,t] +
                    b().K[k].vl * l[i,j,k,t]
                )
                for i in b().K[k].cover, j in b().K[k].cover
            ) +
            sum(
                b().K[k].fd * u[i,k,t]
                for i in b().K[k].cover
            ) +
            sum(
                b().K[k].fp * z[i,k,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                (u[i,k,t] - v[i,k,t]) * duals.λ[i,t]
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
        sum(u[i,k,t] - v[i,k,t] for i in b().K[k].cover) == 0 #all pickup delivered
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        u[i,k,t] <= b().K[k].Q * y[i,k,t] #u-y corr
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        v[i,k,t] <= b().K[k].Q * z[i,k,t] #v-z corr
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(x[j,i,k,t] for j in b().K[k].cover) == y[i,k,t] + z[i,k,t] #traverse in
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(x[i,j,k,t] for j in b().K[k].cover) == y[i,k,t] + z[i,k,t] #traverse out
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(l[j,i,k,t] for j in b().K[k].cover) -
        sum(l[i,j,k,t] for j in b().K[k].cover) == u[i,k,t] - v[i,k,t] #load balance
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, j = b().K[k].cover, t = b().T],
        l[i,j,k,t] <= b().K[k].Q * x[i,j,k,t] #l-x corr
    )

    for k in keys(b().K), t in b().T
        for i in keys(b().V)
            if !(i in b().K[k].cover)
                @constraint(sp, z[i,k,t] == 0)
                @constraint(sp, y[i,k,t] == 0)
            end
        end
    end

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

    return sp
end

function getCols(sp::Model)
    u = value.(sp.obj_dict[:u])
    v = value.(sp.obj_dict[:v])
    l = value.(sp.obj_dict[:l])
    y = value.(sp.obj_dict[:y])
    z = value.(sp.obj_dict[:z])
    x = value.(sp.obj_dict[:x])

    return col(u,v,l,y,z,x)
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
                    println("price: $(objective_value(sp))")
                end

                if isapprox(objective_value(sp),0,atol = 1e-8) || objective_value(sp) > 0
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
                    push!(n.columns,getCols(sp)) #action
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
