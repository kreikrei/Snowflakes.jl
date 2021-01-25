# =========================================================================
#    ENVIRONMENT SETTINGS
# =========================================================================

const default_optimizer = Ref{Any}(nothing)

"""
    set_default_optimizer!(O)
Sets the default optimizer. For example,
    using GLPK
    set_default_optimizer!(GLPK.Optimizer)
"""

set_default_optimizer!(O) = default_optimizer[] = O

"""
    get_default_optimizer()
Gets the default optimizer, which is set by `set_default_optimizer`.
"""

get_default_optimizer() = default_optimizer[]

"""
    reset_default_optimizer()
Resets the default optimizer.
"""

reset_default_optimizer() = default_optimizer[] = nothing
