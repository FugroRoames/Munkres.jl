module Munkres

export munkres

immutable Location
    row::Int
    col::Int
end

Location() = Location(-1,-1)

valid(loc::Location) = loc.row > 0 && loc.col > 0

const StarMark=1
const PrimeMark=2

type ModifiedCost
    m::Matrix{Float64}
    row_offsets::Vector{Float64}
    column_offsets::Vector{Float64}
end

function ModifiedCost{T<:Real}(cost_matrix::Matrix{T})
    ModifiedCost(cost_matrix,
                 zeros(eltype(cost_matrix), size(cost_matrix,1)),
                 zeros(eltype(cost_matrix), size(cost_matrix,2)))
end

@inline Base.size(cost::ModifiedCost, args...) = size(cost.m, args...)

@inline function iszero(cost::ModifiedCost, i, j) 
    @inbounds return cost.m[i,j] == cost.row_offsets[i] + cost.column_offsets[j]
end

@inline function Base.getindex(cost::ModifiedCost, i, j) 
    @inbounds return cost.m[i,j] - cost.row_offsets[i] - cost.column_offsets[j]
end

function find_zeros(cost::ModifiedCost)
    zero_locations = [Int[] for i=1:size(cost,1)]
    for j = 1:size(cost,2)
        for i = 1:size(cost,1)
            @inbounds if iszero(cost, i, j)
                push!(zero_locations[i], j)
            end
        end
    end
    zero_locations
end


"""
Julia implementation of the Munkres (Hungarian) Algorithm for optimal assignment, i.e.
given an n x m cost matrix for assigning n workers to m jobs, what is the optimal
allocation of tasks to workers such that one one task is assigned to one worker.
Based on c# implementation at
"Munkres' Assignment Algorithm, Modified for Rectangular Matrices",
http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
"""
function munkres(cost_matrix)
    n,m = size(cost_matrix)
    row_cover = zeros(Bool,n)
    column_cover = zeros(Bool,m)
    mask_array = zeros(Int8,n,m)
    path_start = Location()    
    cost = ModifiedCost(cost_matrix)

    step_one!(cost)
    step_two!(cost, mask_array, row_cover, column_cover)

    zero_locations = find_zeros(cost)   

    step = 3

    while true
        # println("step = ",step)
        if step == 3
            step = step_three!(mask_array, column_cover)
        elseif step == 4
            step, path_start = step_four!(cost, mask_array, row_cover, column_cover, zero_locations)
        elseif step == 5
            step = step_five!(mask_array, n, m, row_cover, column_cover, path_start)
            path_start = Location()
        elseif step == 6
            step = step_six!(cost,row_cover,column_cover, zero_locations)
            #zero_locations = find_zeros(cost)            
        elseif step == 7
            break
        else
            error("Step has gone bonkers")
        end
        #awesome_print_debug(cost, mask_array, row_cover, column_cover, step, path_start)
    end
    return [findfirst(mask_array[i,:] .== StarMark) for i=1:n]
end

function awesome_print_debug(cost, mask_array, row_cover, column_cover, step, location)
    println("cost matrix = ")
    @show cost.m
    println("mask array = ")
    println(mask_array)
    println("row cover = ")
    println(row_cover)
    println("column cover = ")
    println(column_cover)
    println("augmenting path start = $(location.row) $(location.col)")
    readline(STDIN)
end

function step_one!(cost)
    #remove row minimum from cost matrix
    cost.row_offsets = vec(minimum(cost.m,2))
end

function step_two!(cost, mask_array, row_cover, column_cover)
    #find a zero in cost matrix and star the zero if there is no other zero in row or column
    #(i.e. set mask_array to 1)
    for i=1:size(cost,1)
        for j=1:size(cost,2)
            if iszero(cost, i, j) && !row_cover[i] && !column_cover[j]
                mask_array[i,j] = StarMark
                row_cover[i] = true
                column_cover[j] = true
            end
        end
    end
    row_cover[:] = false
    column_cover[:] = false
end

function step_three!(mask_array, column_cover)
    #cover each column with a starred zero. If k (minimum(n,m)) columns are covered the algorithm is done
    n,m = size(mask_array)
    column_count = 0
    step = 3
    for i=1:n
        for j=1:m
            if mask_array[i,j] == StarMark
                column_cover[j] = true
            end
        end
    end
    column_count = sum(column_cover)
    if (column_count >= m) || (column_count >= n)
       step = 7
    else
       step = 4
    end    
    return step
end

function step_four!(cost, mask_array, row_cover, column_cover, zero_locations)
    # Find indices of all zeros in the uncovered rows.  Doing this as a prework
    # step is faster than doing it in the loop below, probably because the
    # memory traversal order is more cache friendly.
    #zero_locations = find_zeros(cost)

    for row=1:size(cost,1)
        if row_cover[row]
            continue
        end
        # Find column of first zero in cost matrix which is uncovered
        zero_column = 0
        for j = zero_locations[row]
            if !column_cover[j]
                zero_column = j
            end
        end
        if zero_column == 0
            continue
        end

        mask_array[row, zero_column] = PrimeMark
        starcolumn = 0
        for j = 1:size(mask_array,2)
            @inbounds if mask_array[row,j] == StarMark
                starcolumn = j
                break
            end
        end
        if starcolumn == 0
            return 5, Location(row,zero_column)
        end

        row_cover[row] = true
        column_cover[starcolumn] = false
    end

    return 6, Location()
end


function step_five!(mask_array, n, m, row_cover, column_cover, path_start)
    @assert valid(path_start)

    done = false
    row = -1
    column = -1

    path = Array(Location,0)
    push!(path, path_start)

    while ~done
        row = find_star_in_column(mask_array, row, path[end].col, n)
        if row > -1
            push!(path, Location(row, path[end].col))
        else
            done = true
        end
        if ~done
            column = find_prime_in_row(mask_array, path[end].row, column, m)
            push!(path, Location(path[end].row, column))
        end
    end
    augment_path!(path, mask_array)
    row_cover[:] = false
    column_cover[:] = false
    erase_primes!(mask_array)
    return 3
end


function step_six!(cost,row_cover,column_cover, zero_locations)
    min_value, min_locations = find_smallest_uncovered(cost,row_cover,column_cover)
        
    #min locations in the uncovered rows become the new zeros    
    for i=1:length(min_locations)               
        if !isempty(min_locations[i])
            for j=1:length(min_locations[i])
                push!(zero_locations[i],min_locations[i][j])
            end
        end
    end
    
    cost.row_offsets[row_cover] -= min_value
    cost.column_offsets[!column_cover] += min_value
    
    #need to deal with any zeros going away in covered columns and rows
    for i = 1:length(zero_locations)
        if !isempty(zero_locations[i]) && row_cover[i]
            for j=1:length(zero_locations[i])
                if column_cover[zero_locations[i][j]]                                        
                    zero_location = zero_locations[i]
                    zero_locations[i] = zero_location[zero_location .!= zero_locations[i][j]]                    
                    break
                end
            end
        end
    end    
    
    return 4
end


"""
Utility functions for step 5
"""

function find_star_in_column(mask_array, row, active_column, n)
    row = -1
    for i = 1:n
        if mask_array[i,active_column] == StarMark
            row = i
        end
    end
    return row
end

function find_prime_in_row(mask_array, active_row, column, m)
    column = -1
    for j = 1:m
        if mask_array[active_row,j] == PrimeMark
            column = j
        end
    end
    return column
end

function augment_path!(path, mask_array)
    for p = 1:length(path)
        mask_array[path[p].row, path[p].col] = (mask_array[path[p].row, path[p].col] == StarMark) ? 0 : StarMark
    end
end

function erase_primes!(mask_array)
    for j=1:size(mask_array,2)
        for i=1:size(mask_array,1)
            if mask_array[i,j] == PrimeMark
                mask_array[i,j] = 0
            end
        end
    end
end

"""
Utility function for step six
"""

function find_smallest_uncovered(cost, row_cover, column_cover)
    min_value = Inf
    uncovered_row_inds = find(!row_cover);
    
    for j = 1:size(cost,2)
        if !column_cover[j]
            for i = uncovered_row_inds
                @inbounds c = cost[i,j]
                min_value = ifelse(min_value > c, c, min_value)
            end
        end
    end

    min_locations = [Int[] for i=1:size(cost,1)]
    for j = 1:size(cost,2)
        if !column_cover[j]
            for i = uncovered_row_inds
                @inbounds c = cost[i,j]
                if c == min_value
                    push!(min_locations[i],j)
                end
            end
        end
    end

    return min_value, min_locations
end


end # module
