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
        θ[r=keys(R), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T] >= 0
    ) #θ definition
    @variable(mp,
        I[i=keys(d().V), t=vcat(first(d().T)-1,d().T)]
    ) #I definition
    @variable(mp,
        0 <= slack[i=keys(d().V), t=d().T] <= n.stab.slLim[i,t]
    ) #slack
    @variable(mp,
        0 <= surp[i=keys(d().V), t=d().T] <= n.stab.suLim[i,t]
    ) #surplus

    @objective(
        mp, Min,
        sum(
            θ[r,k,f,t] * (
                sum(
                    d().dist[i,j] * (
                        d().K[k].vx * R[r].x[i,j,k,f,t] +
                        d().K[k].vl * R[r].l[i,j,k,f,t]
                    ) for i in d().K[k].cover, j in d().K[k].cover
                ) +
                sum(d().K[k].fd * R[r].u[i,k,f,t] for i in d().K[k].cover) +
                sum(d().K[k].fp * R[r].z[i,k,f,t] for i in d().K[k].cover)
            ) for r in keys(R), k in keys(d().K), f in collect(1:d().K[k].freq), t in d().T
        ) +
        sum(
            d().V[i].h * I[i,t]
            for i in keys(d().V), t in d().T
        ) +
        sum(
            n.stab.slCoeff * slack[i,t]
            for i in keys(d().V), t in d().T
        ) -
        sum(
            n.stab.suCoeff * surp[i,t]
            for i in keys(d().V), t in d().T
        )
    )

    @constraint(
        mp, λ[i=keys(d().V),t=d().T],
        I[i,t-1] + sum(
            R[r].q[i,k,f,t] * θ[r,k,f,t]
            for r in keys(R), k in keys(d().K), f in collect(1:d().K[k].freq)
        ) + slack[i,t] - surp[i,t] == d().d[i,t] + I[i,t]
    )

    @constraint(
        mp, γ[i=keys(d().V),k=keys(d().K),t=d().T],
        sum(
            R[r].z[i,k,f,t] * θ[r,k,f,t]
            for r in keys(R), f in collect(1:d().K[k].freq)
        ) <= d().K[k].freq
    )

    @constraint(
        mp, δ[k=keys(d().K),f=collect(1:d().K[k].freq),t=d().T],
        sum(
            θ[r,k,f,t] for r in keys(R)
        ) <= 1
    )

    @constraint(
        mp, [i=keys(d().V),t=d().T],
        d().V[i].MIN <= I[i,t] <= d().V[i].MAX
    )

    @constraint(
        mp, [i=keys(d().V)],
        I[i,first(d().T)-1] == d().V[i].START
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

    @variable(sp, q[i=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T])
    @variable(sp, u[i=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T] >= 0)
    @variable(sp, v[i=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T] >= 0)
    @variable(sp, l[i=keys(d().V), j=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T] >= 0)
    @variable(sp, p[i=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T], Bin)
    @variable(sp, y[i=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T], Bin)
    @variable(sp, z[i=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T], Bin)
    @variable(sp, x[i=keys(d().V), j=keys(d().V), k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T], Bin)

    @objective(
        sp, Min,
        sum(
            sum(
                d().dist[i,j] * (
                    d().K[k].vx * x[i,j,k,f,t] +
                    d().K[k].vl * l[i,j,k,f,t]
                )
                for i in d().K[k].cover, j in d().K[k].cover
            ) +
            sum(d().K[k].fd * u[i,k,f,t] for i in d().K[k].cover) +
            sum(d().K[k].fp * z[i,k,f,t] for i in d().K[k].cover)
            for k in keys(d().K), f in collect(1:d().K[k].freq), t in d().T
        ) -
        sum(
            q[i,k,f,t] * duals.λ[i,t]
            for k in keys(d().K), i in d().K[k].cover, f in collect(1:d().K[k].freq), t in d().T
        ) -
        sum(
            z[i,k,f,t] * duals.γ[i,k,t]
            for k in keys(d().K), i in d().K[k].cover, f in collect(1:d().K[k].freq), t in d().T
        )
    )

    @constraint(
        sp, [k=keys(d().K), i=d().K[k].cover, f=collect(1:d().K[k].freq), t=d().T],
        q[i,k,f,t] == u[i,k,f,t] - v[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(d().K), i=d().K[k].cover, f=collect(1:d().K[k].freq), t=d().T],
        p[i,k,f,t] == y[i,k,f,t] + z[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(d().K), i=d().K[k].cover, f=collect(1:d().K[k].freq), t=d().T],
        sum(l[j,i,k,f,t] for j in keys(d().V)) -
            sum(l[i,j,k,f,t] for j in keys(d().V)) ==
            q[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(d().K), i=d().K[k].cover, f=collect(1:d().K[k].freq), t=d().T],
        sum(x[j,i,k,f,t] for j in keys(d().V)) +
            sum(x[i,j,k,f,t] for j in keys(d().V)) ==
            2 * p[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(d().K), i=d().K[k].cover, f=collect(1:d().K[k].freq), t=d().T],
        u[i,k,f,t] <= d().K[k].Q * y[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(d().K), i=d().K[k].cover, f=collect(1:d().K[k].freq), t=d().T],
        v[i,k,f,t] <= d().K[k].Q * z[i,k,f,t]
    )

    @constraint(
        sp, [k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T],
        sum(q[i,k,f,t] for i in d().K[k].cover) == 0
    )

    @constraint(
        sp, [k=keys(d().K), f=collect(1:d().K[k].freq), t=d().T],
        sum(z[s,k,f,t] for s in d().K[k].cover) <= 1
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
