module ContinuedFractions
# It'st most appropriate to define a Julia iterable object for this task
# Julia doesn't have Python'st yield, the closest to it is produce/consume calls with Julia tasks
# but for various reasons they don't work out for this task
# This solution works with two integers, a Julia rational or a real

export ContinuedFraction, convergents, intermediate_convergents

mutable struct ContinuedFraction{T<:Integer}
    n1::T # numerator or real
    n2::T # denominator or 1 if real
    t1::T # generated coefficient
end

# Constructors for all possible input types
ContinuedFraction(n1::T, n2::T) where T<:Integer = ContinuedFraction(n1, n2, zero(T))
ContinuedFraction(n::Rational) = ContinuedFraction(numerator(n), denominator(n))
ContinuedFraction(n::AbstractFloat) = ContinuedFraction(Rational(n))

# Outdated Julia 0.6 code
# # Methods to make our object iterable
# Base.start(::ContinuedFraction) = nothing
# # Returns true if we've prepared the continued fraction
# Base.done(cf::ContinuedFraction, st) = cf.n2 == 0
# # Generates the next coefficient
# function Base.next(cf::ContinuedFraction, st)
#     cf.n1, (cf.t1, cf.n2) = cf.n2, divrem(cf.n1, cf.n2)
#     return cf.t1, nothing
# end

function Base.iterate(cf::ContinuedFraction, state=nothing)

    if cf.n2 == 0
        return nothing
    end
    cf.n1, (cf.t1, cf.n2) = cf.n2, divrem(cf.n1, cf.n2)
    return cf.t1, nothing

end


# Tell Julia that this object always returns ints (all coeffs are integers)
Base.eltype(cf::ContinuedFraction) = Integer

# Overload the default collect function so that we can collect the first maxiter coeffs of infinite continued fractions
# array slicing doesn't work as Julia crashes before the slicing due to our infinitely long array
function Base.collect(itr::ContinuedFraction, maxiter::Integer = 100)
    r = Array{eltype(itr)}(undef,maxiter)
    i = 1
    for v in itr
        r[i] = v
        i += 1
        if i > maxiter break end
    end
    return r[1:i-1]
end

function convergents(cf::ContinuedFraction,n::Integer)
    as = collect(cf,n)
    f = Array{Rational}(undef,length(as))
    f[1] = as[1]//1
    if length(as)>=2
        f[2] = (as[2]*as[1]+1)//as[2]
        for j=3:length(as)
            f[j] = (as[j]*numerator(f[j-1]) + numerator(f[j-2]))//(as[j]*denominator(f[j-1]) + denominator(f[j-2]))
        end
    end
    return f
end

function intermediate_convergents(cf::ContinuedFraction,n::Integer)
    as = collect(cf,n)
    n = length(as)
    f = Array{Rational}(undef,2)
    f[1] = as[1]//1
    if n>=2
        f[2] = (as[2]*as[1]+1)//as[2]
        h1 = numerator(f[2])
        h0 = numerator(f[1])
        k1 = denominator(f[2])
        k0 = denominator(f[1])
        for j=3:n
            f = [f ; [(m*h1+h0)//(m*k1+k0) for m=1:as[j]]]
            h0 = h1
            h1 = numerator(f[end])
            k0 = k1
            k1 = denominator(f[end])
        end
    end
    return f
end

end
