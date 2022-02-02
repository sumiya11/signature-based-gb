
include("module.jl")


#-----------------------------------------------------------------------------
# Regular reductions

# try to reduce f with g in a regular way and return result
function regular_reduction_step(vF, f::ModuleElement, g::ModuleElement)
    evalf = vF(f)
    evalg = vF(g)

    @info "reducing $(f) with $(g)" evalf evalg

    for t in terms(evalf)
        # if divides
        success, u = divides(t, leading_term(evalg))
        if success
            u = divexact(t, leading_term(evalg))
            println("$(leading_term(evalg)) | $t")

            # if reduction is regular (and does not change signature)
            # not really correct!
            if signature(f) > signature(u*g)
                newf = f - u*g
                println("regular! $(f) - $(u)*$(g) --> $newf")
                return true, newf
            end
        end
    end
    return false, f
end

# // same as above, but in terms of set G //
# try to reduce f with G in a regular way and return result
function regular_reduction_step(vF, f::ModuleElement, G)
    evalf = vF(f)
    for g in G
        success, newf = regular_reduction_step(vF, f, g)
        if success
            return true, newf
        end
    end
    return false, f
end

# regular normal form of f w.r.t. G
function regular_normal_form(vF, f::ModuleElement, G)
    @info "computing reg normalform of $f w.r.t $G.."

    success = true
    newf = copy(f)
    while success
        success, newf = regular_reduction_step(vF, newf, G)
        @info "reduction $success"
    end
    newf
end

#-----------------------------------------------------------------------------
# Singular reductions

# if f is singurarly top reducible by G
function issingularlytopreducible(vF, f::ModuleElement, G)
    leadevalf = leading_monomial(vF(f))
    for g in G
        evalg = vF(g)
        success, u = divides(leadevalf, leading_monomial(evalg))
        if success # TODO: need to check signature also!!! # UPD: done
            if g*u ≃ f
                return true
            end
        end
    end
    return false
end

#-----------------------------------------------------------------------------
# S-polynomial

# returns u, v,  such that
# u*a - v*b is spoly(a, b)
function mults(vF, a::ModuleElement, b::ModuleElement)
    evala = vF(a)
    evalb = vF(b)

    u = lcm(leading_monomial(evala), leading_monomial(evalb))
    u = divexact(u, leading_term(evala))

    v = lcm(leading_monomial(evala), leading_monomial(evalb))
    v = divexact(v, leading_term(evalb))

    u, v
end

# S-polynomial of a and b
function spoly(vF, a::ModuleElement, b::ModuleElement)
    @info "computing spoly for $a ~ $(vF(a)) and $b ~ $(vF(b))"
    u, v = mults(vF, a, b)
    @info "multipliers are $u and $v: $(u)*$(a) - $(v)*$(b) = $(u*a - v*b)"
    u*a - v*b
end

#-----------------------------------------------------------------------------
# Groebner basis things

# return the i'th basis of element of P^m
function basis_element(F, i)
    R = parent(first(F))
    singlecomponent(one(R), i, length(F))
end

function signature_groebner_basis(F)
    F = map(f -> map_coefficients(c -> c // leading_coefficient(f), f), F)

    m  = length(F)
    vF = x -> v(F, x)

    # fill G
    # create an array of m elements, where i'th element is basis_element(F, i)
    G = [basis_element(F, i) for i in 1:m]

    # fill P
    P = []
    for fi in G
        for fj in G
            u, v = mults(vF, fi, fj)
            if signature(u*fi) >ₚₒₜ signature(v*fj)
                push!(P, spoly(vF, fi, fj))
                # (fi, fj), (fj, fi) ∈ P
            end
        end
    end

    @warn "generated initial G and P:"
    println("F = $F")
    println("G = $G")
    println("P = $P")

    while !isempty(P)
        sigmin = minimum(signature, P)
        println(sigmin)
        index = findfirst(p -> signature(p) ≃ sigmin, P)
        f = P[index] # f is our Spoly
        deleteat!(P, index)
        @warn "selected $f"

        f′ = regular_normal_form(vF, f, G)
        @warn "computed normal form $f′" issyzygy(vF, f′) issingularlytopreducible(vF, f′, G)

        if !issyzygy(vF, f′) && !issingularlytopreducible(vF, f′, G)
            # update P
            for fj in G
                # note that here we take f, and not f′
                u, v = mults(vF, f, fj)
                if !(signature(u*f) ≃ signature(v*fj))
                    # spoly(vF, f, fj) -->
                    push!(P, spoly(vF, f, fj))
                end
            end

            # update G
            f′normalzed = f′ * inv(leading_coefficient(vF(f′)))
            push!(G, f′normalzed)
        end

        @warn "updated G and P"
        println("G = $G")
        println("P = $P")
    end

    G
end

#-----------------------------------------------------------------------------

R, (x,y,z) = PolynomialRing(QQ, ["x","y", "z"], ordering=:degrevlex)

F = [x*y*z + 5, x + y + z, x*y + x*z + y*z]
vF = t -> v(F, t)

# evaluation
G = signature_groebner_basis(F)

println("############")
println(G)
println([vF(g) for g in G])
