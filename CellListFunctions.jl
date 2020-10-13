#
#  CellListFunctions.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 12/10/2020.
#
#

# Functions from https://discourse.julialang.org/t/cell-list-algorithm-is-slower-than-double-for-loop/36621/2

module CellListFunctions

    using LinearAlgebra

    const Intx2 = Tuple{Int64, Int64}
    const Intx3 = Tuple{Int64, Int64, Int64}

    @inline function findPairs!(N, pos, s, neighbourCells)

        pairsList = Intx2[] # Array of tuples storing neighbour pairs
        boundaryList = Intx3[]

        # Allocate all particles in matrix pos to grid points
        cellLists = gridAllocate(pos,N,s)

        for key in keys(cellLists)
            sameCellPairs!(cellLists[key], pairsList, pos, s)
            loopNeighbour!(key, neighbourCells, cellLists, pairsList, pos, s)

            # for jj=1:3
        	# 	if key[jj]==1
        	# 		for kk in cellLists[key]
        	# 			push!(boundaryList,(kk,jj,1))
        	# 		end
        	# 	end
        	# 	if key[jj]==N_grid
        	# 		for kk in cellLists[key]
        	# 			push!(boundaryList,(kk,jj,N_grid))
        	# 		end
        	# 	end
        	# end
        end

        return pairsList#,boundaryList

    end

    # Allocate all particles in matrix pos to grid points. Return Dictionary cellLists mapping (x,y,z) indices to list of particles.
    @inline function gridAllocate(pos::Array{Float64,2}, N::Int64, s::Float64)

        cellLists = Dict{Intx3,Vector{Int64}}()

        for i = 1:N
            indexTuple = (ceil.(Int64,pos[i,:]./s)...,)
            if indexTuple ∉ keys(cellLists)
                cellLists[indexTuple] = [i]
            else
                push!(cellLists[indexTuple],i)
            end
        end

        return cellLists
    end

    # Find list of neighbouring cells around given cell index
    @inline function getForwardCell!(neighbourCells::Vector{Intx3}, index::Intx3)
        neighbourCells[1]  = (index[1]+1, index[2]  , index[3]  )
        neighbourCells[2]  = (index[1]+1, index[2]  , index[3]-1)
        neighbourCells[3]  = (index[1]+1, index[2]  , index[3]+1)
        neighbourCells[4]  = (index[1]+1, index[2]+1, index[3]  )
        neighbourCells[5]  = (index[1]+1, index[2]+1, index[3]-1)
        neighbourCells[6]  = (index[1]+1, index[2]+1, index[3]+1)
        neighbourCells[7]  = (index[1]+1, index[2]-1, index[3]  )
        neighbourCells[8]  = (index[1]+1, index[2]-1, index[3]-1)
        neighbourCells[9]  = (index[1]+1, index[2]-1, index[3]+1)
        neighbourCells[10] = (index[1]  , index[2]  , index[3]+1)
        neighbourCells[11] = (index[1]  , index[2]+1, index[3]  )
        neighbourCells[12] = (index[1]  , index[2]+1, index[3]+1)
        neighbourCells[13] = (index[1]  , index[2]+1, index[3]-1)

        return neighbourCells
    end

    @inline function sameCellPairs!(cellList::Array{Int64}, pairsList::Vector{Intx2}, pos::Matrix{Float64}, s::Float64)
        len = length(cellList)
        if len >= 2
            for i=1:len-1
                for j=i+1:len
                    addOrNot!(cellList[i], cellList[j], pos, s, pairsList)
                end
            end
        end
        return pairsList
    end

    #
    @inline function loopNeighbour!(key::Intx3, neighbourCells::Vector{Intx3}, cellLists::Dict{Intx3,Vector{Int64}}, pairsList::Vector{Intx2}, pos::Matrix{Float64},s ::Float64)
        getForwardCell!(neighbourCells, key)
        for i = 1:4
            if neighbourCells[i] ∈ keys(cellLists)
                findPair!(cellLists[key], cellLists[neighbourCells[i]], pairsList, pos, s)
            end
        end
        return pairsList
    end

    # Function to loop over all particles in two given grid points and then call function to test whether each particle pair is in interaction range
    @inline function findPair!(cellList1::Array{Int64}, cellList2::Array{Int64}, pairsList::Vector{Intx2}, pos::Matrix{Float64}, s::Float64)
        for P1 in cellList1
            for P2 in cellList2
                addOrNot!(P1, P2, pos, s, pairsList)
            end
        end
        return pairsList
    end

    # Test whether two points P1 and P2 are within separation limit s; if so, add tuple (P1,P2) to vector of pairs pairsList
    @inline function addOrNot!(P1::Int64, P2::Int64, pos::Matrix{Float64}, s::Float64, pairsList::Vector{Intx2})
        dx = pos[P1,:].-pos[P2,:]
        distSq = dot(dx,dx)
        if distSq <= s^2
            push!(pairsList, (P1,P2)) # Add this pair to array of neighbour pairs if within range
        end
        return pairsList
    end

export findPairs!

end
