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

    @variable(mp, θ[r=keys(R),k=keys(base().K),t=base().T] >= 0) #θ definition
    @variable(mp, I[i=keys(base().V),t=vcat(first(base().T)-1,base().T)]) #I definition
    @variable(mp, 0 <= slack[i=keys(base().V),t=base().T] <= n.stab.slLim[i,t]) #slack
    @variable(mp, 0 <= surp[i=keys(base().V),t=base().T] <= n.stab.suLim[i,t]) #surplus

    @objective(
        mp, Min,
        #sum(
        #    sum(
        #        base().dist[i,j] * (base().K[k].vx * R[r].x[i,j,k,t] + base().K[k].vl * R[r].l[i,j,k,t]) for i in keys(base().V), j in keys(base().V)) + sum(base().K[k].fd * R[r].u[i,k,t] for i in base().K[k].cover) + sum(base().K[k].fp * R[r].z[i,k,t] for i in base().K[k].loadp) for r in keys(R), k in keys(base().K), t in base().T
        #) +
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
            for i in keys(base().V, t in base().T)
        )
    )

    @constraint(
        mp, λ[i=keys(base().V),t=base().T],
        I[i,t-1] + sum(
            R[r].q[i,k,t] * θ[r,k,t] for r in keys(R), k in keys(base().K)
        ) == base().d[i,t] + I[i,t]
    )

    @constraint(
        mp, γ[i=keys(base().V),k=keys(base().K),t=base().T],
        sum(
            R[r].z[i,k,t] * θ[r,k,t] for r in keys(R)
        ) <= base().K[k].freq
    )

    @constraint(
        mp, δ[k=keys(base().K),t=base().T],
        sum(
            θ[r,k,t] for r in keys(R)
        ) <= base().K[k].freq
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

    return duals
end

function sub(n::node,duals::dval;silent::Bool)

    return sp
end

function buildSub(n::node,duals::dval;silent::Bool)

    return sp
end

function getCols(sp::Model)

    return cols
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
