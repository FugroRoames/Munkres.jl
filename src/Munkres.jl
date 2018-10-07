module Munkres

export munkres

struct Location
    row::Int
    col::Int
end

Location() = Location(-1,-1)

valid(loc::Location) = loc.row > 0 && loc.col > 0

const StarMark=1
const PrimeMark=2

mutable struct ModifiedCost
    m::Matrix{Float64}
    row_offsets::Vector{Float64}
    column_offsets::Vector{Float64}
end

function ModifiedCost(cost_matrix::Matrix{T}) where T<:Real
    ModifiedCost(cost_matrix,
                 zeros(eltype(cost_matrix), size(cost_matrix,1)),
                 zeros(eltype(cost_matrix), size(cost_matrix,2)))
end

Base.size(cost::ModifiedCost, args...) = size(cost.m, args...)

Base.eltype(cost::ModifiedCost) = eltype(cost.m)

function iszero(cost::ModifiedCost, i, j)
    @inbounds return cost.m[i,j] == cost.row_offsets[i] + cost.column_offsets[j]
end

function Base.getindex(cost::ModifiedCost, i, j)
    @inbounds return cost.m[i,j] - cost.row_offsets[i] - cost.column_offsets[j]
end


"""
Munkres (Hungarian) Algorithm for optimal assignment, i.e.  given an NxM cost
matrix for assigning N workers to M jobs, what is the optimal allocation of
tasks to workers such that one one task is assigned to one worker.

Input

  * cost_matrix - an NxM real valued matrix of the cost of assigning job m to worker N

Output

  * p - A vector of length N of assignments, where j = p[i] assigns job j to worker i.

"""
function munkres(cost_matrix)
    # Inspired by Bob Pilgrim's C# tutorial, "Munkres' Assignment Algorithm,
    # Modified for Rectangular Matrices",
    # http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
    # and by and Yi Cao's surprisingly fast matlab version,
    # http://mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm

    any(isnan, cost_matrix) && error("cost matrix cannot have NaNs")
    n,m = size(cost_matrix)
    flipped = false
    if n > m
        #always use more jobs than workers (current implmentation doesn't work in other case, so just transpose)
        cost_matrix = Array(cost_matrix')
        flipped = true
        n,m = size(cost_matrix)
    end
    row_cover = zeros(Bool,n)
    column_cover = zeros(Bool,m)
    mask_array = zeros(Int8,n,m)
    path_start = Location()
    cost = ModifiedCost(cost_matrix)

    zero_locations = step_one!(cost)
    step_two!(cost, mask_array, row_cover, column_cover)

    step = 3

    while true
        #println("step = ",step)
        if step == 3
            step = step_three!(mask_array, column_cover)
        elseif step == 4
            step, path_start = step_four!(cost, mask_array, row_cover, column_cover, zero_locations)
        elseif step == 5
            step = step_five!(mask_array, row_cover, column_cover, path_start)
            path_start = Location()
        elseif step == 6
            step = step_six!(cost,row_cover,column_cover, zero_locations)
        elseif step == 7
            break
        else
            error("Step has gone bonkers")
        end
    end

    if flipped
        return [findfirst(mask_array[:,i] .== StarMark) for i=1:size(mask_array,2)]
    end
    return [findfirst(mask_array[i,:] .== StarMark) for i=1:size(mask_array,1)]
end


function step_one!(cost)
    #remove row minimum from cost matrix and find locations of all zeros
    cost.row_offsets = vec(minimum(cost.m, dims=2))
    zero_locations = [findall(j->iszero(cost,i,j), 1:size(cost,2)) for i=1:size(cost,1)]
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
    row_cover[:] .= false
    column_cover[:] .= false
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


function step_five!(mask_array, row_cover, column_cover, path_start)
    #Find a series of alternating starred and primed zeros in alternating columns and rows
    #until it terminates on a primed zero with no star in it's column. Unstar all the starred
    #zeros, star the primed zeros, erase the primes and uncover all lines.
    @assert valid(path_start)

    done = false
    row = -1
    column = -1

    path = Vector{Location}()
    push!(path, path_start)

    while ~done
        row = find_star_in_column(mask_array, row, path[end].col)
        if row > -1
            push!(path, Location(row, path[end].col))
        else
            done = true
        end
        if ~done
            column = find_prime_in_row(mask_array, path[end].row, column)
            push!(path, Location(path[end].row, column))
        end
    end
    augment_path!(path, mask_array)
    row_cover[:] .= false
    column_cover[:] .= false
    erase_primes!(mask_array)
    return 3
end


function step_six!(cost,row_cover,column_cover, zero_locations)
    #Subtract the minimum value in the uncovered part of the cost matrix from all uncovered
    #columns and add it to the covered rows. Note that for efficiency this does not explicitly
    #occur, we just keep track of where new zeros appear and where zeros disappear and update
    #the row and column minima.
    min_value, min_locations = find_smallest_uncovered(cost,row_cover,column_cover)

    #min locations in the uncovered rows become the new zeros
    for i=1:length(min_locations)
        push!(zero_locations[min_locations[i][1]],min_locations[i][2])
    end

    cost.row_offsets[row_cover] .-= min_value
    cost.column_offsets[map(!, column_cover)] .+= min_value

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


function find_star_in_column(mask_array, row, active_column)
    #find the row in which a star occurs for the active column
    row = -1
    for i = 1:size(mask_array,1)
        if mask_array[i,active_column] == StarMark
            row = i
        end
    end
    return row
end


function find_prime_in_row(mask_array, active_row, column)
    #find the column in which a prime occurs for the active row
    column = -1
    for j = 1:size(mask_array,2)
        if mask_array[active_row,j] == PrimeMark
            column = j
        end
    end
    return column
end


function augment_path!(path, mask_array)
    #go along the path swapping stars for zeros
    for p = 1:length(path)
        mask_array[path[p].row, path[p].col] = (mask_array[path[p].row, path[p].col] == StarMark) ? 0 : StarMark
    end
end


function erase_primes!(mask_array)
    #remove all of the primes from the mask_array
    for j=1:size(mask_array,2)
        for i=1:size(mask_array,1)
            if mask_array[i,j] == PrimeMark
                mask_array[i,j] = 0
            end
        end
    end
end

function find_smallest_uncovered(cost, row_cover, column_cover)
    #find the locations and value of the minimum of the cost matrix in the uncovered rows and columns
    min_value = typemax(eltype(cost))
    uncovered_row_inds = findall(map(!, row_cover))
    uncovered_col_inds = findall(map(!, column_cover))
    min_locations = Tuple{Int, Int}[]
    for j in uncovered_col_inds, i in uncovered_row_inds
        @inbounds c = cost[i,j]
        if c < min_value
            #have found a new minimum, so throw away all the previously discovered values and start again
            min_value = c
            empty!(min_locations)
            push!(min_locations,(i,j))
        elseif  c == min_value
            push!(min_locations,(i,j))
        end
    end
    return min_value, min_locations
end

end # module
