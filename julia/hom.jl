function homogenize(F)
    base = parent(F[1])
    s = copy(symbols(base))
    push!(s, :h)
    R, xsu = PolynomialRing(base_ring(base), s, ordering=ordering(base))
    xs = xsu[1:end-1]
    u = xsu[end]

    homF = Vector{eltype(F)}(undef, 0)
    for f in F
        tdeg = total_degree(leading_monomial(f))
        f = evaluate(f, xs)
        homf = zero(R)
        for t in terms(f)
            homf += t * u^(tdeg - total_degree(t))
        end
        push!(homF, homf)
    end

    homF
end

function dehomogenize(homF)
    R = parent(homF[1])
    xs = collect(gens(R))
    xs[end] = one(R)

    F = Vector{eltype(homF)}(undef, 0)
    for homf in homF
        f = evaluate(homf, xs)
        push!(F, f)
    end

    F
end
