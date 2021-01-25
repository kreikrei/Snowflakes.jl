# ==============================================================================
#    COLUMN GENERATION MECHANISMS
# ==============================================================================

function master(n::node;silent::Bool)

    return mp
end

function buildMaster(n::node;silent::Bool)

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
