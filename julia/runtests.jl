
include("groebner.jl")
using Groebner

###############################################################################

fs = Groebner.rootn(3)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
eval_sign_gb = v(fs, sign_gb)

println(eval_sign_gb)
@assert Groebner.isgroebner(eval_sign_gb)

###############################################################################

fs = Groebner.change_ordering(Groebner.noonn(2), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
eval_sign_gb = v(fs, sign_gb)

println(eval_sign_gb)
@assert Groebner.isgroebner(eval_sign_gb)
