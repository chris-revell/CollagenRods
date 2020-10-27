#
#  ShortestDistance.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module ShortestDistance

using LinearAlgebra
include("./ShortestPointToRod.jl")
using .ShortestPointToRod

@inline function shortestDistance!(r,Ω,rᵢⱼ,i,j,L)
    # Algorithm from Vega, Lago 1994
    rᵢⱼ .= r[j,:] .- r[i,:]
    λᵢ = (rᵢⱼ⋅Ω[i,:] - (Ω[i,:]⋅Ω[j,:])*(rᵢⱼ⋅Ω[j,:]))/(1.0-(Ω[i,:]⋅Ω[j,:])^2)
    μⱼ = ((Ω[i,:]⋅Ω[j,:])*(rᵢⱼ⋅Ω[i,:]) - rᵢⱼ⋅Ω[j,:])/(1.0-(Ω[i,:]⋅Ω[j,:])^2)

    if -L/2.0 <= λᵢ <= L/2.0 && -L/2.0 <= μⱼ <= L/2.0
        # Shortest path is within length of rods and not from an end of either rod
        λ = λᵢ
        μ = μⱼ
    else
        # Shortest path involves a rod end
        θ = atan(λᵢ,μⱼ) + π/4.0
        if 0 < θ <= π/2.0
            # Region 1
            μ = L/2.0
            Pⱼ = r[j,:].+μ.*Ω[j,:]
            t = shortestPointToRod(Pⱼ,r[i,:],Ω[i,:],L)
            λ = (t-0.5)*L
        elseif π/2.0 < θ <= π
            λ = L/2.0
            Pᵢ = r[i,:] .+ λ.*Ω[i,:]
            t = shortestPointToRod(Pᵢ,r[j,:],Ω[j,:],L)
            μ = (t-0.5)*L
        elseif π < θ <= 3.0*π/2.0
            μ = -L/2.0
            Pⱼ = r[j,:].+μ.*Ω[j,:]
            t = shortestPointToRod(Pⱼ,r[i,:],Ω[i,:],L)
            λ = (t-0.5)*L
        else
            λ = -L/2.0
            Pᵢ = r[i,:] .+ λ.*Ω[i,:]
            t = shortestPointToRod(Pᵢ,r[j,:],Ω[j,:],L)
            μ = (t-0.5)*L
        end
    end

    rᵢⱼ .+= (μ.*Ω[j,:] .- λ.*Ω[i,:])
    return (μ,λ)
end

export shortestDistance!

end
