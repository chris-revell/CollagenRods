#
#  CellListFunctions.jl
#  collagen-rods
#
#  Created by Christopher Revell on 12/10/2020.
#
#
# Functions to allocate particle positions to cell lists in a grid, and hence produce a list of neighbour pairs.
# Derived from https://discourse.julialang.org/t/cell-list-algorithm-is-slower-than-double-for-loop/36621/2

module CellListFunctions

    using LinearAlgebra
    using Dictionaries

    const Intx2 = Tuple{Int64, Int64}
    const Intx3 = Tuple{Int64, Int64, Int64}

    @inline @views function findPairs!(N, r, interactionThresh, neighbourCells)

        pairsList = Intx2[] # Array of tuples storing neighbour pairs
        #boundaryList = Intx3[]

        # Allocate all particles in matrix r to grid points
        cellLists = gridAllocate(r,N,interactionThresh)

        for key in keys(cellLists)
            sameCellPairs!(cellLists[key], pairsList, r, interactionThresh)
            loopNeighbour!(key, neighbourCells, cellLists, pairsList, r, interactionThresh)

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

    # Allocate all particles in matrix r to grid points. Return Dictionary cellLists mapping (x,y,z) indices to list of particles.
    @inline @views function gridAllocate(r, N, interactionThresh)

        cellLists = Dictionary{Intx3,Vector{Int64}}()

        for i = 1:N
            indexTuple = (ceil.(Int64,r[i,:]./interactionThresh)...,)
            if indexTuple ∉ keys(cellLists)
                insert!(cellLists,indexTuple,[i])
            else
                push!(cellLists[indexTuple],i)
            end
        end

        return cellLists
    end

    # Find list of neighbouring cells around given cell index
    @inline @views function getForwardCell!(neighbourCells, index)
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

    @inline @views function sameCellPairs!(cellList, pairsList, r, interactionThresh)
        len = length(cellList)
        if len >= 2
            for i=1:len-1
                for j=i+1:len
                    addOrNot!(cellList[i], cellList[j], r, interactionThresh, pairsList)
                end
            end
        end
        return pairsList
    end

    #
    @inline @views function loopNeighbour!(key, neighbourCells, cellLists, pairsList, r,interactionThresh)
        getForwardCell!(neighbourCells, key)
        for neighbourIndex in neighbourCells
            if neighbourIndex ∈ keys(cellLists)
                findPair!(cellLists[key], cellLists[neighbourIndex], pairsList, r, interactionThresh)
            end
        end
        return pairsList
    end

    # Function to loop over all particles in two given grid points and then call function to test whether each particle pair is in interaction range
    @inline @views function findPair!(cellList1, cellList2, pairsList, r, interactionThresh)
        for P1 in cellList1
            for P2 in cellList2
                addOrNot!(P1, P2, r, interactionThresh, pairsList)
            end
        end
        return pairsList
    end

    # Test whether two points P1 and P2 are within separation limit interactionThresh; if so, add tuple (P1,P2) to vector of pairs pairsList
    @inline @views function addOrNot!(P1, P2, r, interactionThresh, pairsList)
        dx = r[P1,:].-r[P2,:]
        distSq = dx⋅dx
        if distSq <= interactionThresh^2
            push!(pairsList, (P1,P2)) # Add this pair to array of neighbour pairs if within range
        end
        return pairsList
    end

export findPairs!

end
