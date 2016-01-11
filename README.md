# Munkres

Julia implementation of the
[Hungarian algorithm](https://en.wikipedia.org/wiki/Hungarian_algorithm)
for the optimal assignment problem: Given N workers and M jobs, find a minimal
cost assignment of one job to each worker.

[![Build Status](https://travis-ci.org/FugroRoames/Munkres.jl.svg?branch=master)](https://travis-ci.org/FugroRoames/Munkres.jl)

## Examples

A synthetic example with a simple solution.

```julia
# Each worker-job combination has a random cost
cost = rand(4,4)
# However, each worker can do a certain job with zero cost
best_jobs = [3,4,1,2]
for (i,j) in enumerate(best_jobs); cost[i,j] = 0; end

# Compute optimal assignment given the cost
computed_best_jobs = munkres(cost)

@assert best_jobs == computed_best_jobs
```


Example output:

```julia
julia> cost = rand(4,4)
4x4 Array{Float64,2}:
 0.455632  0.0972016  0.807122  0.806295
 0.933324  0.280094   0.261727  0.235289
 0.53323   0.408037   0.935853  0.62243
 0.08281   0.147279   0.649129  0.910296

julia> best_jobs = [3,4,1,2]
4-element Array{Int64,1}:
 3
 4
 1
 2

julia> for (i,j) in enumerate(best_jobs); cost[i,j] = 0; end

julia> computed_best_jobs = munkres(cost)
4-element Array{Int64,1}:
 3
 4
 1
 2

```
