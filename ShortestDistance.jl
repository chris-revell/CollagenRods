#
#  ShortestDistance.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module ShortestDistance

using LinearAlgebra

@inline function shortestRodToRod!(r,Ω,rᵢⱼ,i,j,L)
    # Algorithm from Vega, Lago 1994
    rᵢⱼ .= view(r,j,:) .- view(r,i,:)
    λᵢ = (rᵢⱼ⋅view(Ω,i,:) - (view(Ω,i,:)⋅view(Ω,j,:))*(rᵢⱼ⋅view(Ω,j,:)))/(1.0-(view(Ω,i,:)⋅view(Ω,j,:))^2)
    μⱼ = ((view(Ω,i,:)⋅view(Ω,j,:))*(rᵢⱼ⋅view(Ω,i,:)) - rᵢⱼ⋅view(Ω,j,:))/(1.0-(view(Ω,i,:)⋅view(Ω,j,:))^2)

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
            Pⱼ = view(r,j,:).+μ.*view(Ω,j,:)
            t = shortestPointToRod(Pⱼ,r[i,:],Ω[i,:],L)
            λ = (t-0.5)*L
        elseif π/2.0 < θ <= π
            λ = L/2.0
            Pᵢ = view(r,i,:) .+ λ.*view(Ω,i,:)
            t = shortestPointToRod(Pᵢ,r[j,:],Ω[j,:],L)
            μ = (t-0.5)*L
        elseif π < θ <= 3.0*π/2.0
            μ = -L/2.0
            Pⱼ = view(r,j,:).+μ.*view(Ω,j,:)
            t = shortestPointToRod(Pⱼ,r[i,:],Ω[i,:],L)
            λ = (t-0.5)*L
        else
            λ = -L/2.0
            Pᵢ = view(r,i,:) .+ λ.*view(Ω,i,:)
            t = shortestPointToRod(Pᵢ,r[j,:],Ω[j,:],L)
            μ = (t-0.5)*L
        end
    end

    rᵢⱼ .+= (μ.*view(Ω,j,:) .- λ.*view(Ω,i,:))
    return (μ,λ)
end

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
export shortestRodToRod!

end
