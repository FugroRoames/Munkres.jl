using Munkres
using Base.Test


"""
Solve assignment problem via brute force exploration over all permutations
"""
function brute_force_optima(cost_matrix)
    (n,m) = size(cost_matrix)
    if n != m
        error("Non-square cost matrix, I haven't gotten around to doing that yet!")
    end
    initial_permuation = 1:n
    best_permuation = 1:n
    min_cost = Inf
    for p in permutations(initial_permuation)
        current_cost = 0.0
        for i = 1:n
            current_cost += cost_matrix[initial_permuation[i],p[i]]
        end
        if current_cost <= min_cost
            min_cost = current_cost
            best_permuation = p
        end
    end
    return best_permuation
end


# Small matrix designed to hit all steps
tst = [1 2 3;
       2 4 6;
       3 6 9]
@test munkres(tst) == brute_force_optima(tst)


for i=1:10
    tst = rand(8,8)
    @test munkres(tst) == brute_force_optima(tst)
end

