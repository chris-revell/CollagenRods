#
#  cellListFunctions.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 03/08/2020.
#
#

# Functions from https://discourse.julialang.org/t/cell-list-algorithm-is-slower-than-double-for-loop/36621/2

module CellListFunctions

    using LinearAlgebra

    const Intx2 = Tuple{Int64, Int64}
    const Intx3 = Tuple{Int64, Int64, Int64}

    @inline function find_pairs(N_particles::Int64, pos::Array{Float64,2}, s::Float64, N_grid::Int64, neighbour_cells::Vector{Intx3})

        pairs_list = Intx2[] # Array of tuples storing neighbour pairs
        boundary_list = Intx3[]

        # Allocate all particles in matrix pos to grid points
        cell_lists = grid_allocate(pos,N_particles,s)

        for key in keys(cell_lists)
            same_cell_pairs!(cell_lists[key], pairs_list, pos, s)
            loop_neighbor!(key, neighbour_cells, cell_lists, pairs_list, pos, s)

            for jj=1:3
        		if key[jj]==1
        			for kk in cell_lists[key]
        				push!(boundary_list,(kk,jj,1))
        			end
        		end
        		if key[jj]==N_grid
        			for kk in cell_lists[key]
        				push!(boundary_list,(kk,jj,N_grid))
        			end
        		end
        	end

        end

        return pairs_list#,boundary_list

    end

    # Allocate all particles in matrix pos to grid points. Return Dictionary cell_lists mapping (x,y,z) indices to list of particles.
    @inline function grid_allocate(pos::Array{Float64,2}, N_particles::Int64, s::Float64)

        cell_lists = Dict{Intx3,Vector{Int64}}()

        for i = 1:N_particles
            index_tuple = (ceil.(Int64,pos[i,:]./s)...,)
            if index_tuple ∉ keys(cell_lists)
                cell_lists[index_tuple] = [i]
            else
                push!(cell_lists[index_tuple],i)
            end
        end

        return cell_lists
    end

    # Find list of neighbouring cells around given cell index
    @inline function get_forward_cell!(neighbour_cells::Vector{Intx3}, index::Intx3)
        neighbour_cells[1]  = (index[1]+1, index[2]  , index[3]  )
        neighbour_cells[2]  = (index[1]+1, index[2]  , index[3]-1)
        neighbour_cells[3]  = (index[1]+1, index[2]  , index[3]+1)
        neighbour_cells[4]  = (index[1]+1, index[2]+1, index[3]  )
        neighbour_cells[5]  = (index[1]+1, index[2]+1, index[3]-1)
        neighbour_cells[6]  = (index[1]+1, index[2]+1, index[3]+1)
        neighbour_cells[7]  = (index[1]+1, index[2]-1, index[3]  )
        neighbour_cells[8]  = (index[1]+1, index[2]-1, index[3]-1)
        neighbour_cells[9]  = (index[1]+1, index[2]-1, index[3]+1)
        neighbour_cells[10] = (index[1]  , index[2]  , index[3]+1)
        neighbour_cells[11] = (index[1]  , index[2]+1, index[3]  )
        neighbour_cells[12] = (index[1]  , index[2]+1, index[3]+1)
        neighbour_cells[13] = (index[1]  , index[2]+1, index[3]-1)

        return neighbour_cells
    end

    @inline function same_cell_pairs!(cell_list::Array{Int64}, pairs_list::Vector{Intx2}, pos::Matrix{Float64}, s::Float64)
        len = length(cell_list)
        if len >= 2
            for i=1:len-1
                for j=i+1:len
                    add_or_not!(cell_list[i], cell_list[j], pos, s, pairs_list)
                end
            end
        end
        return pairs_list
    end

    #
    @inline function loop_neighbor!(key::Intx3, neighbour_cells::Vector{Intx3}, cell_lists::Dict{Intx3,Vector{Int64}}, pairs_list::Vector{Intx2}, pos::Matrix{Float64},s ::Float64)
        get_forward_cell!(neighbour_cells, key)
        for i = 1:4
            if neighbour_cells[i] ∈ keys(cell_lists)
                find_pair!(cell_lists[key], cell_lists[neighbour_cells[i]], pairs_list, pos, s)
            end
        end
        return pairs_list
    end

    # Function to loop over all particles in two given grid points and then call function to test whether each particle pair is in interaction range
    @inline function find_pair!(cell_list1::Array{Int64}, cell_list2::Array{Int64}, pairs_list::Vector{Intx2}, pos::Matrix{Float64}, s::Float64)
        for P1 in cell_list1
            for P2 in cell_list2
                add_or_not!(P1, P2, pos, s, pairs_list)
            end
        end
        return pairs_list
    end

    # Test whether two points P1 and P2 are within separation limit s; if so, add tuple (P1,P2) to vector of pairs pairs_list
    @inline function add_or_not!(P1::Int64, P2::Int64, pos::Matrix{Float64}, s::Float64, pairs_list::Vector{Intx2})
        dx = pos[P1,:].-pos[P2,:]
        dist_sq = dot(dx,dx)
        if dist_sq <= s^2
            push!(pairs_list, (P1,P2)) # Add this pair to array of neighbour pairs if within range
        end
        return pairs_list
    end

export find_pairs

end
