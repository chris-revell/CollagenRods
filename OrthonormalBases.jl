#
#  OrthonormalBases.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module OrthonormalBases

using LinearAlgebra

@inline function orthonormalBases!(N,Ω,E)
    for i=1:N
        # Create orthonormal basis vectors around rod axis
        E[i,:,:] .= nullspace(Matrix(Ω[i,:]'))
    end
    return nothing
end

export orthonormalBases!

end
