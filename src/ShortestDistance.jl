#
#  ShortestDistance.jl
#  CollagenRods
#
#  Created by Christopher Revell on 07/10/2020.
#
#
# Functions to find the shortest distances between two rods and between a point and a rod. 

module ShortestDistance

using LinearAlgebra

@inline @views function shortestRodToRod!(r,Ω,rᵢⱼ,i,j,L,dummyVectors)
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
            dummyVectors[1,:] = r[j,:].+μ.*Ω[j,:]
            t = shortestPointToRod(dummyVectors[1,:],dummyVectors[2,:],dummyVectors[3,:],r[i,:],Ω[i,:],L)
            λ = (t-0.5)*L
        elseif π/2.0 < θ <= π
            λ = L/2.0
            dummyVectors[1,:] = r[i,:] .+ λ.*Ω[i,:]
            t = shortestPointToRod(dummyVectors[1,:],dummyVectors[2,:],dummyVectors[3,:],r[j,:],Ω[j,:],L)
            μ = (t-0.5)*L
        elseif π < θ <= 3.0*π/2.0
            μ = -L/2.0
            dummyVectors[1,:] = r[j,:].+μ.*Ω[j,:]
            t = shortestPointToRod(dummyVectors[1,:],dummyVectors[2,:],dummyVectors[3,:],r[i,:],Ω[i,:],L)
            λ = (t-0.5)*L
        else
            λ = -L/2.0
            dummyVectors[1,:] = r[i,:] .+ λ.*Ω[i,:]
            t = shortestPointToRod(dummyVectors[1,:],dummyVectors[2,:],dummyVectors[3,:],r[j,:],Ω[j,:],L)
            μ = (t-0.5)*L
        end
    end

    rᵢⱼ .+= (μ.*Ω[j,:] .- λ.*Ω[i,:])
    return (μ,λ)
end

@inline @views function shortestPointToRod(x₀,x₁,x₂,r,Ω,L)
    # Algorithm from https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    # Returns 3 vector location of point on rod with centrepoint r, length L, and orientation Ω that is closest to point x₀

    x₁ .= r .- (L/2.0).*Ω
    x₂ .= r .+ (L/2.0).*Ω

    t = -(x₁.-x₀)⋅Ω/(L.*Ω⋅Ω)

    if t<0
        t = 0
    elseif t>1
        t = 1
    else
        t = t
    end

    return t

end

export shortestPointToRod
export shortestRodToRod!

end
