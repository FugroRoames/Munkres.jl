module Munkres

export munkres


immutable Location
    row::Int
    col::Int
end

Location() = Location(-1,-1)

valid(loc::Location) = loc.row > 0 && loc.col > 0

"""
Julia implementation of the Munkres (Hungarian) Algorithm for optimal assignment, i.e.
given an n x m cost matrix for assigning n workers to m jobs, what is the optimal
allocation of tasks to workers such that one one task is assigned to one worker.
Based on c# implementation at
"Munkres' Assignment Algorithm, Modified for Rectangular Matrices",
http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
"""
function munkres(cost_matrix)
    initial_cost_matrix = cost_matrix
    cost_matrix = copy(cost_matrix)

    n,m = size(cost_matrix)
    row_cover = zeros(Bool,n)
    column_cover = zeros(Bool,m)
    mask_array = zeros(Int8,n,m)
    path_start = Location()

    step_one!(cost_matrix)
    step_two!(cost_matrix, mask_array, row_cover, column_cover)

    step = 3

    while true
        # println("step = ",step)
        if step == 3
            step = step_three!(mask_array, column_cover)
        elseif step == 4
            step, path_start = step_four!(cost_matrix, mask_array, row_cover, column_cover)
        elseif step == 5
            step = step_five!(mask_array, n, m, row_cover, column_cover, path_start)
            path_start = Location()
        elseif step == 6
            step = step_six!(cost_matrix,row_cover,column_cover)
        elseif step == 7
            break
        else
            error("Step has gone bonkers")
        end
        # awesome_print_debug(cost_matrix, mask_array, row_cover, column_cover, step, path_start)
    end
    return [findfirst(mask_array[i,:] .== 1) for i=1:n]
end

function awesome_print_debug(cost_matrix, mask_array, row_cover, column_cover, step, location)
    println("cost matrix = ")
    println(cost_matrix)
    println("mask array = ")
    println(mask_array)
    println("row cover = ")
    println(row_cover)
    println("column cover = ")
    println(column_cover)
    println("augmenting path start = $(location.row) $(location.col)")
    readline(STDIN)
end

function step_one!(cost_matrix)
    #remove row minimum from cost matrix
    row_min = minimum(cost_matrix,2)
    cost_matrix[:] = cost_matrix .- row_min
end

function step_two!(cost_matrix, mask_array, row_cover, column_cover)
    #find a zero in cost matrix and star the zero if there is no other zero in row or column
    #(i.e. set mask_array to 1)
    for i=1:size(cost_matrix,1)
        for j=1:size(cost_matrix,2)
            if cost_matrix[i,j] == 0.0 && !row_cover[i] && !column_cover[j]
                mask_array[i,j] = 1
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
            if mask_array[i,j] == 1
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

function step_four!(cost_matrix, mask_array, row_cover, column_cover)
    while true
        location = find_a_zero(cost_matrix, row_cover, column_cover)
        row = location.row
        column = location.col
        if (row == -1)
            return 6, Location()
        else
            mask_array[row, column] = 2
            if star_in_row(mask_array, row)
                column = find_star_in_row(mask_array, column, row)
                row_cover[row] = true
                column_cover[column] = 0
            else
                return 5, Location(row,column)
            end
        end
    end
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
    mask_array = erase_primes(mask_array, n, m)
    return 3
end

function step_six!(cost_matrix,row_cover,column_cover)
    min_value = find_smallest_uncovered(cost_matrix,row_cover,column_cover)
    for j=1:size(cost_matrix,2)
        for i=1:size(cost_matrix,1)
            if row_cover[i] == true
                cost_matrix[i,j] += min_value
            end
            if column_cover[j] == false
                cost_matrix[i,j] -= min_value
            end
        end
    end
    step = 4
end


"""
Utility functions for step 4
"""

function find_a_zero(cost_matrix, row_cover, column_cover)
    for j = 1:size(cost_matrix,2), i = 1:size(cost_matrix,1)
        @inbounds if (cost_matrix[i,j] == 0.0) && (!row_cover[i]) && (!column_cover[j])
            return Location(i,j)
        end
    end
    return Location(-1,-1)
end

function star_in_row(mask_array, row)
    star = false
    for j = 1:size(mask_array,2)
        if mask_array[row,j] == 1
            star = true
        end
    end
    return star
end

function find_star_in_row(mask_array, column, row)
    column = -1
    for j = 1:size(mask_array,2)
        if mask_array[row,j] == 1
            column = j
        end
    end
    return column
end

"""
Utility functions for step 5
"""

function find_star_in_column(mask_array, row, active_column, n)
    row = -1
    for i = 1:n
        if mask_array[i,active_column] == 1
            row = i
        end
    end
    return row
end

function find_prime_in_row(mask_array, active_row, column, m)
    column = -1
    for j = 1:m
        if mask_array[active_row,j] == 2
            column = j
        end
    end
    return column
end

function augment_path!(path, mask_array)
    for p = 1:length(path)
        mask_array[path[p].row, path[p].col] = (mask_array[path[p].row, path[p].col] == 1) ? 0 : 1
    end
end

function erase_primes(mask_array, n, m)
    for i=1:n
        for j=1:m
            if mask_array[i,j] == 2
                mask_array[i,j] = 0
            end
        end
    end
    return mask_array
end

"""
Utility function for step six
"""

function find_smallest_uncovered(cost_matrix, row_cover, column_cover)
    min_value = Inf
    for j=1:size(cost_matrix,2)
        for i=1:size(cost_matrix,1)
            @inbounds if !row_cover[i] && !column_cover[j]
                current_value = cost_matrix[i,j]
                if min_value > cost_matrix[i,j]
                    min_value = cost_matrix[i,j]
                end
            end
        end
    end
    return min_value
end


end # module
