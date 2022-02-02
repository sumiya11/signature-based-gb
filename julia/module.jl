
# import some useful things

import Base: +, -, *, <
import Base: copy, show
# provides polynomial ring and polynomials interface
using AbstractAlgebra
# provides dot product operation
import LinearAlgebra: dot

#-----------------------------------------------------------------------------

# our module element will store one thing: a tuple of polynomials
struct ModuleElement
    mtuple
end

# some magic to allow copying ModuleElement
function Base.copy(m::ModuleElement)
    return ModuleElement(deepcopy(m.mtuple))
end

# some magic to print ModuleElement nicely
Base.show(io::IO, ::MIME"text/plain", m::ModuleElement) = println(io, "($(join(map(string, m.mtuple), ", ")))")
Base.show(io::IO, m::ModuleElement) = print(io, "($(join(map(string, m.mtuple), ", ")))")

#-----------------------------------------------------------------------------
# Orderings stuff

# compare two module elements w.r.t. pot strategy,
# a nice thing about julia is that we can use some latex capabilities in naming variables
function <ₚₒₜ(a::ModuleElement, b::ModuleElement)
    # compare indices of first nonzero
    ia = findfirst(!iszero, a.mtuple)
    ib = findfirst(!iszero, b.mtuple)
    if ia != ib
        return ia > ib
    end

    # compare leading monomials if indices are equal,
    # `leading_monomial` comes from `AbstractAlgebra`
    leadcompa = leading_monomial(a.mtuple[ia])
    leadcompb = leading_monomial(b.mtuple[ib])
    leadcompa < leadcompb
end

>ₚₒₜ(a::ModuleElement, b::ModuleElement) = b <ₚₒₜ a

<(a::ModuleElement, b::ModuleElement) = a <ₚₒₜ b
Base.isless(a::ModuleElement, b::ModuleElement) = a < b

#-----------------------------------------------------------------------------
# Some stuff

# creates (0, 0, f,...,0, 0) of length n where f is at position i
function singlecomponent(f, i::Int, n::Int)
    mtuple = zeros(parent(f), n)
    mtuple[i] = f
    ModuleElement(mtuple)
end

# extend leading_monomial definition to our module element
function AbstractAlgebra.leading_monomial(a::ModuleElement)
    ia = findfirst(!iszero, a.mtuple)
    # if all components are zeros
    ia == nothing && return singlecomponent(zero(a.mtuple[1]), 1, length(a.mtuple))

    f = leading_monomial(a.mtuple[ia])
    singlecomponent(f, ia, length(a.mtuple))
end

function AbstractAlgebra.leading_term(a::ModuleElement)
    ia = findfirst(!iszero, a.mtuple)
    ia == nothing && return singlecomponent(zero(a.mtuple[1]), 1, length(a.mtuple))

    f = leading_term(a.mtuple[ia])
    singlecomponent(f, ia, length(a.mtuple))
end

#-----------------------------------------------------------------------------
# Arithmetic operations

# returns t*m where t is something, and m is a module element
*(t, m::ModuleElement) = m * t
*(m::ModuleElement, t) = ModuleElement(m.mtuple .* t)

# group operations
-(a::ModuleElement) = ModuleElement(- a.mtuple)
+(a::ModuleElement, b::ModuleElement) = ModuleElement(a.mtuple .+ b.mtuple)
-(a::ModuleElement, b::ModuleElement) = a + (-b)

function (≃)(a::ModuleElement, b::ModuleElement)
    leading_monomial(a).mtuple == leading_monomial(b).mtuple
end

#-----------------------------------------------------------------------------
# Signatures and evaluation

function signature(a::ModuleElement)
    leading_term(a)
end

dot(F, m::ModuleElement) = dot(m, F)
dot(m::ModuleElement, F) = sum(m.mtuple .* F)

# evaluation hom
v(F, m::ModuleElement) = dot(F, m)
v(F, M) = [v(F, m) for m in M]

issyzygy(vF, h::ModuleElement) = iszero(vF(h))
