#
#  ShortestDistance.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module ShortestDistance

using LinearAlgebra

@inline function shortestDistance!(r,Ω,r₁₂,i,j)
    # Algorithm from Vega, Lago 1994
    r₁₂ .= r[i,:].-r[j,:]
    λᵢ = (r₁₂⋅Ω[i,:] - (Ω[i,:]⋅Ω[j,:])*(r₁₂⋅Ω[j,:]))/(1.0-(Ω[i,:]⋅Ω[j,:])^2)
    μⱼ = ((Ω[i,:]⋅Ω[j,:])*(r₁₂⋅Ω[i,:]) - r₁₂⋅Ω[j,:])/(1.0-(Ω[i,:]⋅Ω[j,:])^2)
    r₁₂ .+= μⱼ.*Ω[j,:] .- λᵢ.*Ω[i,:]
    return (μⱼ,λᵢ)
end

export shortestDistance!

end
