# =========================================================================
#    ENVIRONMENT SETTINGS
# =========================================================================

const default_optimizer = Ref{Any}(nothing)
const slack_coeff = Ref{Any}(nothing)
const surp_coeff = Ref{Any}(nothing)
const verbosity = Ref{Any}(nothing)

set_default_optimizer!(O) = default_optimizer[] = O
get_default_optimizer() = default_optimizer[]
reset_default_optimizer() = default_optimizer[] = nothing

set_slack_coeff!(O) = slack_coeff[] = O
sl_C() = slack_coeff[]

set_surp_coeff!(O) = surp_coeff[] = O
su_C() = surp_coeff[]

set_silent!(O) = verbosity[] = O
silent() = verbosity[]
