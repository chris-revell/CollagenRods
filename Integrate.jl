#
#  Integrate.jl
#  collagen-rods
#
#  Created by Christopher Revell on 10/12/2020.
#
#
# Function to integrate position and orientation of all rods

module Integrate

using LinearAlgebra

@inline @views function integrate!(r,Ω,F,τ,ξr,ξΩ,Δt)

    Ω .+= τ[:,:,1].*Δt .+ ξΩ.*sqrt(Δt)
    normalize!.(eachslice(Ω,dims=1))
    r .+= F[:,:,1].*Δt .+ ξr.*sqrt(Δt)

    return (t + Δt)
end

export integrate!

end
