module TravellingSalesman

using JuMP, HiGHS, Distances

struct Problem
    jumpmodel
    locations
    n::Int
end

# TODO: This should be a generalized module.
# TODO: Create separate functions for "model cycle" and "model path"
function build_tsp_model_for_season(circuits)::Problem
    locations = circuits |> eachrow |> collect
    d = [Distances.haversine((a["lng"], a["lat"]), (b["lng"], b["lat"])) for a ∈ locations, b ∈ locations]

    # Add location that has free travel to serve as our "cut" point to break the cycle into a path.
    d = vcat(hcat(d,
                  zeros(Int, length(locations))),
             zeros(Int, length(locations)+1)')

    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "output_flag", false)
    n = size(d, 1)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @objective(model, Min, sum(d .* x)/2)
    @constraint(model, [i ∈ 1:n], sum(x[i, :]) == 2)
    @constraint(model, [i ∈ 1:n], x[i, i] == 0)
    Problem(model, locations, n)
end

selected_edges(x::Matrix{Float64}, n) = Tuple{Int,Int}[(i, j) for i ∈ 1:n, j ∈ 1:n if x[i, j] > 0.5]

function subtour(edges::Vector{Tuple{Int,Int}}, n)
    shortest_subtour, unvisited = collect(1:n), Set(collect(1:n))
    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors = [j for (i, j) ∈ edges if i == current && j ∈ unvisited]
        end
        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
        end
    end
    shortest_subtour
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

function solve!(prob::Problem)
    optimize!(prob.jumpmodel)
    global cycle = subtour(prob.jumpmodel[:x])
    while 1 < length(cycle) < prob.n
        S = [(i, j) for (i, j) ∈ Iterators.product(cycle, cycle) if i < j]
        @constraint(prob.jumpmodel,
                    sum(prob.jumpmodel[:x][i, j] for (i, j) ∈ S) <= length(cycle) - 1,
                    )
        optimize!(prob.jumpmodel)
        global cycle = subtour(prob.jumpmodel[:x])
    end

    startend = []
    for (i, j) ∈ selected_edges(value.(prob.jumpmodel[:x]), prob.n)
        if i == prob.n
            push!(startend, j)
        elseif j == prob.n
            push!(startend, i)
        end
    end
    unique!(startend)

    grid = value.(prob.jumpmodel[:x])

    orderedindexes = [first(startend)]
    while length(orderedindexes) < prob.n - 1
        row = grid[last(orderedindexes), 1:prob.n-1]
        for (i, v) ∈ row |> enumerate
            if v < 0.5 || i ∈ orderedindexes
                continue
            end

            push!(orderedindexes, i)
            break
        end
    end

    orderedindexes .|> oidx -> prob.locations[oidx]
end

end # module TravellingSalesman
