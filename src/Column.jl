# =========================================================================
#    COLUMN GENERATION MECHANISMS
# =========================================================================

function Q(key,R::Dict)
    if isa(key,β)
        q = Vector{NamedTuple}()

        for r in keys(R), k in keys(b().K), t in b().T
            if key.t == "≳"
                if getproperty(R[r],key.q)[key.i,k,t] >= key.v
                    push!(q,(r=r,k=k,t=t))
                end
            elseif key.t == "<"
                if getproperty(R[r],key.q)[key.i,k,t] < key.v
                    push!(q,(r=r,k=k,t=t))
                end
            end
        end

        return q
    elseif isa(key,Vector{β})
        q = Vector{Vector{NamedTuple}}()

        for b in key
            push!(q,Q(b,R))
        end

        if !isempty(q)
            #ngembaliin semua rkt
            return reduce(intersect,q)
        else
            #ngembaliin vector kosong
            return q
        end
    end
end

function f(B,R,θ)
    if isempty(B)
        return sum(θ[r,k,t] - floor(θ[r,k,t])
            for r in keys(R), k in keys(b().K), t in b().T
        )
    else
        if !isempty(Q(B,R))
            return sum(θ[q.r,q.k,q.t] - floor(θ[q.r,q.k,q.t])
                for q in Q(B,R)
            )
        else
            return 0
        end
    end
end

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

    ≲ = filter(b -> last(b).t == "≲",B)
    ≳ = filter(b -> last(b).t == "≳",B)

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
    @constraint(mp, μ[b = keys(≲)],
        sum(θ[q.r,q.k,q.t] for q in Q(B[b].B,R)) <= B[b].κ
    )

    @constraint(mp, ν[b = keys(≳)],
        sum(θ[q.r,q.k,q.t] for q in Q(B[b].B,R)) >= B[b].κ
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

    ≲ = filter(b -> last(b).t == "≲",B)
    ≳ = filter(b -> last(b).t == "≳",B)

    # ================================
    #    MODEL CONSTRUCTION
    # ================================
    q = col(
        @variable(sp, u[keys(b().V), keys(b().K), b().T] >= 0, Int),
        @variable(sp, v[keys(b().V), keys(b().K), b().T] >= 0, Int),
        @variable(sp,
            l[collect(keys(b().V)), collect(keys(b().V)), collect(keys(b().K)), b().T] >=
            0, Int
        ),
        @variable(sp, y[keys(b().V), keys(b().K), b().T] >= 0, Int),
        @variable(sp, 0 <= z[i = keys(b().V), k = keys(b().K), t = b().T] <=
            b().K[k].freq, Int
        ),
        @variable(sp,
            x[collect(keys(b().V)), collect(keys(b().V)), collect(keys(b().K)), b().T] >=
            0, Int
        )
    )

    @variable(sp, g[keys(≲), keys(b().K), b().T], Bin)
    @variable(sp, h[keys(≳), keys(b().K), b().T], Bin)

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
            sum(
                g[j,k,t] * duals.μ[j]
                for j in keys(≲)
            )
            for k in keys(b().K), t in b().T
        ) -
        sum(
            sum(
                h[j,k,t] * duals.ν[j]
                for j in keys(≳)
            )
            for k in keys(b().K), t in b().T
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
        n = [i for i in keys(b().V) if !(i in b().K[k].cover)] #sets not in cover

        if !isempty(n)
            @constraint(sp, z[n,k,t] == 0)
            @constraint(sp, y[n,k,t] == 0)
        end
    end

    # ================================
    #    BOUND GENERATOR
    # ================================
    imax()

    for j in keys(≲), k in keys(b().K), t in b().T
        η = @variable(sp, [B[j].B], Bin)

        @constraint(sp, -g[j,k,t] >= 1 - sum(1 - η[e] for e in B[j].B))
        for e in B[j].B
            if e.t == "≳"
                @constraint(sp, (getproperty(qmax(),e.q)[e.i] - e.v + 1) * η[e] >=
                    getproperty(q,e.q)[e.i,k,t] - e.v + 1
                )
            elseif e.t == "<"
                @constraint(sp, e.v * η[e] >= e.v - getproperty(q,e.q)[e.i,k,t])
            end
        end
    end

    for j in keys(≳), k in keys(b().K), t in b().T
        η = @variable(sp, [B[j].B], Bin)

        @constraint(sp, [e = B[j].B], h[j,k,t] <= η[e])
        for e in B[j].B
            if e.t == "≳"
                @constraint(sp, e.v * η[e] <= getproperty(q,e.q)[e.i,k,t])
            elseif e.t == "<"
                @constraint(sp,(getproperty(qmax(),e.q)[e.i] - e.v + 1) * η[e] <=
                    (getproperty(qmax(),e.q)[e.i] - getproperty(q,e.q)[e.i,k,t])
                )
            end
        end
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
    s = (
        sum(value.(mp.obj_dict[:slack])) + sum(value.(mp.obj_dict[:surp]))
    )

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

                if isapprox(objective_value(sp),0,atol = 1e-7) || objective_value(sp) > 0
                    if isapprox(checkStab(mp),0,atol = 1e-7)
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
