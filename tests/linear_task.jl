module p

# include("./make_network.jl")
using JuMP
using SCIP
using Gurobi
using Pipe
using OrderedCollections
using PyCall
import JSON

plt = PyCall.pyimport("matplotlib.pyplot")


function calc_dist(x1, x2)
    return(
        @pipe (x2 - x1).^2 |>
        sum(_) |>
        sqrt(_)
    )
end


function print_soln(x, y, u)
    println("x")
    x_soln = @pipe (value.(x) |> round.(Int, _))
    println(x_soln)
    for i in ALL_NODES, j in ALL_NODES
        if x_soln[i, j] == 1
            println("$i -> $j")
        end
    end

    println("\ny")
    y_soln = @pipe (value.(y) |> round.(Int, _))
    println(y_soln)

    println("\nu")
    u_soln = @pipe value.(u) |> round.(_, digits=2)
    for (ui, ui_val) in zip(u.axes[1], u_soln)
        println("$ui = $(ui_val) t units")
    end
end


"""
Plots a solution tour inclduing start, end, defects, compulsory edges and directions
"""
function plot_soln(x, y, u)

    # plot nodes
    xs = [v[1] for (k, v) in points]
    ys = [v[2] for (k, v) in points]
    point_names = [k for (k, v) in points]

    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 10])
    ax.grid(linewidth=0.2)
    ax.set_title("Tour solution")

    # add node labels 
    for ni in DEFECT_NODES
        xi, yi = points[ni]
        ax.plot(xi, yi, "k.")
    end

    xi, yi = points["start"]
    ax.plot(xi, yi, "r^")

    xi, yi = points["end"]
    ax.plot(xi, yi, "rv")

    for ni in ALL_NODES
        xi, yi = points[ni]
        ax.text(xi, yi, ni)
    end

    # plot tour sequence in order of time progression
    sorted_ks = sortperm(value.(u).data)
    sorted_names = value.(u).axes[1][sorted_ks]
    sorted_us = value.(u).data[sorted_ks]
    for k in 1:NUM_POINTS - 1
        v1 = points[sorted_names[k]]
        v2 = points[sorted_names[k + 1]]
        u1 = sorted_us[k]
        u2 = sorted_us[k + 1]

        # direction vector
        dv = v2 - v1

        # only plot edges with a positive time
        if u1[1] == 0
            continue
        end

        ax.plot([v1[1], v2[1]], [v1[2], v2[2]], linestyle=":", color="b", linewidth=0.7)
        ax.annotate("", xy=(v1 + 0.2*dv), xytext=v1, arrowprops=Dict("arrowstyle"=>"->", "color"=>"g"))
    end

    # plot compulsory edges
    for (n1, n2) in compulsory_edges
        v1 = points[n1]
        v2 = points[n2]
        ax.plot([v1[1], v2[1]], [v1[2], v2[2]], "r-")
    end

    plt.show()
end


"""
Returns a network in standard form (as a set of points and compulsory edges) by
reading a json file created using make_network.jl
"""
function read_json_network(json_file::String)
    loaded_network_data = nothing
    open(json_file, "r") do f
        loaded_network_data = readlines(f)[1] |> JSON.parse
    end

    # data type conversions
    points = OrderedDict{String, Vector{Float64}}(loaded_network_data["points"])
    compulsory_edges = Vector{Vector{String}}(loaded_network_data["compulsory_edges"])

    return points, compulsory_edges
end


# ---------------------------------------------------------------------------

println("loading network into model")

# loading json network data
points, compulsory_edges = read_json_network("C:\\Users\\Elton.Shi\\Desktop\\GitHub\\DownerRoads\\routine_inspection_optimisation_project\\tests\\test_network_5.json")

# compulsory_edges = network.compulsory_edges
# points = network.all_pts

# Sets
NUM_POINTS = length(points)
DEFECT_NAMES = [ki for ki in keys(points) if !(ki in ("start", "end"))]
DEFECT_NODES = OrderedSet(DEFECT_NAMES)
START_NODE = OrderedSet(["start"])
END_NODE = OrderedSet(["end"])
ALL_NODES = union(START_NODE, DEFECT_NODES, END_NODE)

D = JuMP.Containers.DenseAxisArray{Float64}(undef, ALL_NODES, ALL_NODES)

# create distance matrix
for i in ALL_NODES, j in ALL_NODES
    D[i, j] = calc_dist(points[i], points[j])
end


# model setup
# -----------
m = Model(Gurobi.Optimizer)

# decision variables
# ------------------
@variable(m, x[ALL_NODES, ALL_NODES], Bin)
@variable(m, y[ALL_NODES], Bin)
@variable(m, u[ALL_NODES], Int)


# constraints
# -----------

# 1. compulsory edges must be done
for (i, j) in compulsory_edges
    @constraint(m, x[i, j] == 1)

    # prohibit reverse travel
    @constraint(m, x[j, i] == 0)
end

# 2. If node j is active, one of the entering edges attached to j must be active
for j in union(DEFECT_NODES, END_NODE)
    @constraint(m, y[j] <= sum(x[i, j] for i in union(START_NODE, DEFECT_NODES)))
end

# 3. if exiting node j, one of the exiting edges must be active 
for i in union(START_NODE, DEFECT_NODES)
    @constraint(m, y[i] <= sum(x[i, j] for j in union(DEFECT_NODES, END_NODE)))
end

# 4. cannot travel to self
for i in ALL_NODES
    @constraint(m, x[i, i] == 0)
end

# 5. if a node i is selected, only one edge can be active leaving the node
for i in union(START_NODE, DEFECT_NODES)
    @constraint(m, sum(x[i, j] for j in union(DEFECT_NODES, END_NODE)) <= 1)
end

# if node j is selected, only one edge can be active entering the node
for j in union(DEFECT_NODES, END_NODE)
    @constraint(m, sum(x[i, j] for i in union(START_NODE, DEFECT_NODES)) <= 1)
end

# 6. Cannot enter start node, cannot exit end node
for i in union(DEFECT_NODES, END_NODE)
    @constraint(m, x[i, "start"] == 0)
end
for j in union(START_NODE, DEFECT_NODES)
    @constraint(m, x["end", j] == 0)
end

# 7. if x_ij is active, y_i and y_j must be active
for i in ALL_NODES, j in ALL_NODES
    @constraint(m, x[i,j] <= y[i])
    @constraint(m, x[i,j] <= y[j])
end

# 8. MTZ subtour constraints: if x_ij is active, t_j >= t_i + margin
MARGIN_GAP = 1
M = length(ALL_NODES)
for i in ALL_NODES, j in ALL_NODES
    @constraint(m, u[i] <= u[j] - MARGIN_GAP + M*(1 - x[i, j]))
end

# start node potential must take the value of 1
@constraint(m, u["start"] == 1.0)

# if a node is not used, it's node potential must be zero
for i in DEFECT_NODES
    @constraint(m, u[i] <= M*y[i])
end

# all nodes potentials cannot be negative
for i in ALL_NODES
    @constraint(m, u[i] >= 0)
end

# debug: limit max number of nodes
# @constraint(m, sum(y[i] for i in ALL_NODES) <= 12)


# objective
# ---------
@expression(m, tot_nodes_visited, sum(y[i] for i in ALL_NODES))
@expression(m, tot_dist_travelled, sum(D[i, j] * x[i, j] for i in ALL_NODES, j in ALL_NODES))
@expression(m, obj_fn, tot_nodes_visited - 0.1*tot_dist_travelled)
@objective(m, Max, obj_fn)


# solve
# -----
optimize!(m)


# display results
# ---------------
print_soln(x, y, u)
plot_soln(x, y, u)

end