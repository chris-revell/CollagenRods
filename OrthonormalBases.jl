#
#  OrthonormalBases.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module OrthonormalBases

using LinearAlgebra
using StaticArrays
using Base.Threads

@inline @views function orthonormalBases!(N,Ω,E)
    @threads for i=1:N
        # Create orthonormal basis vectors around rod axis
        E[i] = SMatrix{3,2}(nullspace(SMatrix{1,3}(Ω[i]')))
    end
    return nothing
end

export orthonormalBases!

end
