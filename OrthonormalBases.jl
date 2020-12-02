#
#  OrthonormalBases.jl
#  collagen-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#
# Function to calculate perpendicular vectors around each rod given that rod's orientation. 

module OrthonormalBases

using LinearAlgebra
using Base.Threads

@inline @views function orthonormalBases!(N,Ω,E)
    @threads for i=1:N
        # Create orthonormal basis vectors around rod axis
        E[i,:,:] .= nullspace(Matrix((Ω[i,:])'))
    end
    return nothing
end

export orthonormalBases!

end
