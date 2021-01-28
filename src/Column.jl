# =========================================================================
#    COLUMN GENERATION MECHANISMS
# =========================================================================

function master(n::node;silent::Bool)
    mp = buildMaster(n;silent = silent)
    optimize!(mp)

    return mp
end

function buildMaster(n::node;silent::Bool)
    R = Dict(1:length(n.columns) .=> n.columns)
    mp = Model(get_default_optimizer())
    if silent
        set_silent(mp)
    end

    @variable(mp,
        θ[r=keys(R), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T] >= 0
    ) #θ definition
    @variable(mp,
        I[i=keys(base().V), t=vcat(first(base().T)-1,base().T)]
    ) #I definition
    @variable(mp,
        0 <= slack[i=keys(base().V), t=base().T] <= n.stab.slLim[i,t]
    ) #slack
    @variable(mp,
        0 <= surp[i=keys(base().V), t=base().T] <= n.stab.suLim[i,t]
    ) #surplus

    @objective(
        mp, Min,
        sum(
            θ[r,k,f,t] * (
                sum(
                    base().dist[i,j] * (
                        base().K[k].vx * R[r].x[i,j,k,f,t] +
                        base().K[k].vl * R[r].l[i,j,k,f,t]
                    ) for i in base().K[k].cover, j in base().K[k].cover
                ) +
                sum(base().K[k].fd * R[r].u[i,k,f,t] for i in base().K[k].cover) +
                sum(base().K[k].fp * R[r].z[i,k,f,t] for i in base().K[k].cover)
            ) for r in keys(R), k in keys(base().K), f in collect(1:base().K[k].freq), t in base().T
        ) +
        sum(
            base().V[i].h * I[i,t]
            for i in keys(base().V), t in base().T
        ) +
        sum(
            n.stab.slCoeff * slack[i,t]
            for i in keys(base().V), t in base().T
        ) -
        sum(
            n.stab.suCoeff * surp[i,t]
            for i in keys(base().V), t in base().T
        )
    )

    @constraint(
        mp, λ[i=keys(base().V),t=base().T],
        I[i,t-1] + sum(
            R[r].q[i,k,f,t] * θ[r,k,f,t]
            for r in keys(R), k in keys(base().K), f in collect(1:base().K[k].freq)
        ) + slack[i,t] - surp[i,t] == base().d[i,t] + I[i,t]
    )

    @constraint(
        mp, γ[i=keys(base().V),k=keys(base().K),t=base().T],
        sum(
            R[r].z[i,k,f,t] * θ[r,k,f,t]
            for r in keys(R), f in collect(1:base().K[k].freq)
        ) <= base().K[k].freq
    )

    @constraint(
        mp, δ[k=keys(base().K),f=collect(1:base().K[k].freq),t=base().T],
        sum(
            θ[r,k,f,t] for r in keys(R)
        ) <= 1
    )

    @constraint(
        mp, [i=keys(base().V),t=base().T],
        base().V[i].MIN <= I[i,t] <= base().V[i].MAX
    )

    @constraint(
        mp, [i=keys(base().V)],
        I[i,first(base().T)-1] == base().V[i].START
    )

    return mp
end

function getDuals(mp::Model)
    λ = dual.(mp.obj_dict[:λ])
    γ = dual.(mp.obj_dict[:γ])
    δ = dual.(mp.obj_dict[:δ])

    return dval(λ,γ,δ)
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

    @variable(sp, q[i=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T])
    @variable(sp, u[i=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T] >= 0)
    @variable(sp, v[i=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T] >= 0)
    @variable(sp, l[i=keys(base().V), j=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T] >= 0)
    @variable(sp, p[i=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T], Bin)
    @variable(sp, y[i=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T], Bin)
    @variable(sp, z[i=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T], Bin)
    @variable(sp, x[i=keys(base().V), j=keys(base().V), k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T], Bin)

    @objective(
        sp, Min,
        sum(
            sum(
                base().dist[i,j] * (
                    base().K[k].vx * x[i,j,k,f,t] +
                    base().K[k].vl * l[i,j,k,f,t]
                )
                for i in base().K[k].cover, j in base().K[k].cover
            ) +
            sum(base().K[k].fd * u[i,k,f,t] for i in base().K[k].cover) +
            sum(base().K[k].fp * z[i,k,f,t] for i in base().K[k].cover)
            for k in keys(base().K), f in collect(1:base().K[k].freq), t in base().T
        ) -
        sum(
            q[i,k,f,t] * duals.λ[i,t]
            for k in keys(base().K), i in base().K[k].cover, f in collect(1:base().K[k].freq), t in base().T
        ) -
        sum(
            z[i,k,f,t] * duals.γ[i,k,t]
            for k in keys(base().K), i in base().K[k].cover, f in collect(1:base().K[k].freq), t in base().T
        )
    )

    @constraint(
        sp, [k=keys(base().K), i=base().K[k].cover, f=collect(1:base().K[k].freq), t=base().T],
        q[i,k,f,t] == u[i,k,f,t] - v[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(base().K), i=base().K[k].cover, f=collect(1:base().K[k].freq), t=base().T],
        p[i,k,f,t] == y[i,k,f,t] + z[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(base().K), i=base().K[k].cover, f=collect(1:base().K[k].freq), t=base().T],
        sum(l[j,i,k,f,t] for j in keys(base().V)) -
            sum(l[i,j,k,f,t] for j in keys(base().V)) ==
            q[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(base().K), i=base().K[k].cover, f=collect(1:base().K[k].freq), t=base().T],
        sum(x[j,i,k,f,t] for j in keys(base().V)) +
            sum(x[i,j,k,f,t] for j in keys(base().V)) ==
            2 * p[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(base().K), i=base().K[k].cover, f=collect(1:base().K[k].freq), t=base().T],
        u[i,k,f,t] <= base().K[k].Q * y[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(base().K), i=base().K[k].cover, f=collect(1:base().K[k].freq), t=base().T],
        v[i,k,f,t] <= base().K[k].Q * z[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T],
        sum(q[i,k,f,t] for i in base().K[k].cover) == 0
    )

    @constraint(
        sp, [k=keys(base().K), f=collect(1:base().K[k].freq), t=base().T],
        sum(z[s,k,f,t] for s in base().K[k].cover) <= 1
    )

    return sp
end

function getCols(sp::Model)
    p = value.(sp.obj_dict[:p])
    q = value.(sp.obj_dict[:q])
    u = value.(sp.obj_dict[:u])
    y = value.(sp.obj_dict[:y])
    v = value.(sp.obj_dict[:v])
    z = value.(sp.obj_dict[:z])
    l = value.(sp.obj_dict[:l])
    x = value.(sp.obj_dict[:x])

    return col(q,u,v,l,p,y,z,x)
end

function updateStab(stab::stabilizer)

    return stab
end

function checkStab(mp::Model)

    return sum
end

function colGen(n::node;silent::Bool,maxCG::Float64,track::Bool)

    return n
end
