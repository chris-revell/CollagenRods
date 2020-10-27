#
#  ShortestPointToRod.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 26/10/2020.
#
#

module ShortestPointToRod

using LinearAlgebra

@inline function shortestPointToRod(x₀,r,Ω,L)
    # Algorithm from https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    # Returns 3 vector location of point on rod with centrepoint r, length L, and orientation Ω that is closest to point x₀

    x₁ = r .- (L/2.0).*Ω
    x₂ = r .+ (L/2.0).*Ω
    Δx = Ω.*L

    t = -(x₁.-x₀)⋅Δx/(Δx⋅Δx)

    if t<0
        t = 0 #P₂ = r .- (L/2.0).*Ω
    elseif t>1
        t = 1 #P₂ = r .+ (L/2.0).*Ω
    else
        t = t #P₂ = r .- (L/2.0).*Ω .+ t.*Δj
    end

    return t

end

export shortestPointToRod

end
