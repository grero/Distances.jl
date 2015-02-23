# Bhattacharyya distances. Much like for KLDivergence we assume the vectors to
# be compared are probability distributions, frequencies or counts rather than
# vectors of samples. Pre-calc accordingly if you have samples.

type CauchySchwartzDist <: Metric end

# Cauchy-Schwartz coefficient

function cauchyschwartz_coeff{T<:Number}(a::AbstractVector{T}, b::AbstractVector{T})
    n = length(a)
    sqab = zero(T)
    # We must normalize since we cannot assume that the vectors are normalized to probability vectors.
    asum = zero(T)
    bsum = zero(T)

    for i = 1:n
        @inbounds ai = a[i]
        @inbounds bi = b[i]
        sqab += ai * bi
        asum += ai*ai
        bsum += bi*bi
    end

    sqab*sqab / (asum * bsum)
end

cauchyschwartz_coeff{T <: Number}(a::T, b::T) = throw("Cauchy-Schwartz coefficient cannot be calculated for scalars")

# Faster pair- and column-wise versions TBD...


# Cauchy-Schwartz distance
evaluate{T<:Number}(dist::CauchySchwartzDist, a::AbstractVector{T}, b::AbstractVector{T}) = -log(cauchyschwartz_coeff(a, b))
cauchyschwartz(a::AbstractVector, b::AbstractVector) = evaluate(CauchySchwartzDist(), a, b)
evaluate{T <: Number}(dist::CauchySchwartzDist, a::T, b::T) = throw("Cauchy-Schwartz distance cannot be calculated for scalars")
cauchyschwartz{T <: Number}(a::T, b::T) = evaluate(CauchySchwartzDist(), a, b)
