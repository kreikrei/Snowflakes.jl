# =========================================================================
#    COLUMN GENERATION MECHANISMS
# =========================================================================

function master(n::node;silent::Bool)
    mp = buildMaster(n;silent = silent)
    optimize!(mp)

    return mp
end

function buildMaster(n::node;silent::Bool)
    mp = Model(get_default_optimizer())
    if silent
        set_silent(mp)
    end

    R = Dict(1:length(n.columns) .=> n.columns)

    @variable(mp, θ[keys(R), keys(b().K), b().T] >= 0) #column decision
    @variable(mp, I[keys(b().V), vcat(first(b().T)-1, b().T)]) #inventory level
    @variable(mp, 0 <= slack[i = keys(b().V), t = b().T] <= n.stab.slLim[i,t]) #slack var
    @variable(mp, 0 <= surp[i = keys(b().V), t = b().T] <= n.stab.suLim[i,t]) #surp var

    @objective(
        mp, Min,
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
        ) + #pickup delivery costs
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
        I[i,t-1] + sum(R[r].q[i,k,t] * θ[r,k,t] for r in keys(R), k in keys(b().K)) +
        slack[i,t] - surp[i,t] == b().d[i,t] + I[i,t] #inventory balance
    )

    @constraint(mp, δ[i = keys(b().V), k = keys(b().K), t = b().T],
        sum(R[r].z[i,k,t] * θ[r,k,t] for r in keys(R)) <= b().K[k].freq
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
    @constraint(mp, [g = n.bounds],
        sum(
            θ[r,g.idx.k,g.idx.t]
            for r in keys(filter(p -> last(p).z[g.idx.i,g.idx.k,g.idx.t] > 0,R))
        ) == g.val
    )

    return mp
end

function getDuals(mp::Model)
    λ = dual.(mp.obj_dict[:λ])
    δ = dual.(mp.obj_dict[:δ])

    return dval(λ,δ)
end

function sub(n::node,duals::dval;silent::Bool)
    sp = buildSub(n,duals;silent=silent)
    optimize!(sp)

    return sp
end

function buildSub(n::node,duals::dval;silent::Bool)
    sp = Model(get_default_optimizer())
    if silent
        set_silent(sp)
    end

    @variable(sp, q[i = keys(b().V), k = keys(b().K), t = b().T])
    @variable(sp, u[i = keys(b().V), k = keys(b().K), t = b().T] >= 0)
    @variable(sp, v[i = keys(b().V), k = keys(b().K), t = b().T] >= 0)
    @variable(sp,
        l[i = collect(keys(b().V)), j = collect(keys(b().V)),
        k = collect(keys(b().K)), t = b().T] >= 0
    ) #quantity variables

    @variable(sp, p[i = keys(b().V), k = keys(b().K), t = b().T] >= 0, Int)
    @variable(sp, y[i = keys(b().V), k = keys(b().K), t = b().T] >= 0, Int)
    @variable(sp, z[i = keys(b().V), k = keys(b().K), t = b().T] >= 0, Int)
    @variable(sp,
        x[i = collect(keys(b().V)), j = collect(keys(b().V)),
        k = collect(keys(b().K)), t = b().T] >= 0, Int
    ) #0-1 variables

    @objective(
        sp, Min,
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
        ) - #column cost
        sum(
            sum(
                q[i,k,t] * duals.λ[i,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) - #dual part 1
        sum(
            sum(
                z[i,k,t] * duals.δ[i,k,t]
                for i in b().K[k].cover
            )
            for k in keys(b().K), t in b().T
        ) #dual part 2
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        q[i,k,t] == u[i,k,t] - v[i,k,t] #q breakdown
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        p[i,k,t] == y[i,k,t] + z[i,k,t] #p breakdown
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        u[i,k,t] <= b().K[k].Q * y[i,k,t] # u - y correlation
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        v[i,k,t] <= b().K[k].Q * z[i,k,t] # v - z correlation
    )

    @constraint(sp,
        [k = keys(b().K), i = b().K[k].cover, j = b().K[k].cover, t = b().T],
        l[i,j,k,t] <= b().K[k].Q * x[i,j,k,t] # l - x correlation
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        z[i,k,t] <= b().K[k].freq #maximum manifest from each point
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(x[j,i,k,t] for j in b().K[k].cover) == p[i,k,t] #traverse in to i
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(x[i,j,k,t] for j in b().K[k].cover) == p[i,k,t] #traverse out from i
    )

    @constraint(sp, [k = keys(b().K), i = b().K[k].cover, t = b().T],
        sum(l[j,i,k,t] for j in b().K[k].cover) -
        sum(l[i,j,k,t] for j in b().K[k].cover) == q[i,k,t] #vehicle load balance
    )

    @constraint(sp, [k = keys(b().K), t = b().T],
        sum(q[i,k,t] for i in b().K[k].cover) == 0 #all pickup delivered
    )

    return sp
end

function getCols(sp::Model)
    q = value.(sp.obj_dict[:q])
    u = value.(sp.obj_dict[:u])
    v = value.(sp.obj_dict[:v])
    l = value.(sp.obj_dict[:l])
    p = value.(sp.obj_dict[:p])
    y = value.(sp.obj_dict[:y])
    z = value.(sp.obj_dict[:z])
    x = value.(sp.obj_dict[:x])

    return col(q,u,v,l,p,y,z,x)
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

function colGen(n::node;silent::Bool,maxCG::Float64,track::Bool)
    terminate = false
    iter = 0

    while !terminate
        if iter < maxCG
            mp = master(n;silent=silent)
            if has_values(mp) && has_duals(mp)
                if track #print master problem obj
                    println("obj: $(objective_value(mp))")
                end

                duals = getDuals(mp)
                sp = sub(n,duals;silent=silent)
                if track #print subproblem price
                    println("price: $(objective_value(sp))")
                end

                if isapprox(objective_value(sp), 0, atol = 1e-8) || objective_value(sp) > 0
                    if isapprox(checkStab(mp), 0, atol = 1e-8)
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
            push!(n.status,"EVALUATED-TIME OUT") #report
            if track
                println("EVALUATED-TIME OUT")
            end
        end
    end

    if n.status != "NO_SOLUTION"
        println("NODE $(n.self) FINISHED.")
    else
        println("NODE $(n.self) FAILED.")
    end

    return n
end

function matrix(n::node)

    return convert
end
