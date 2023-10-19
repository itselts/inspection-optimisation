module network
using JuMP
using SCIP
using Pipe
import PyPlot
using PyCall
pygui(true)
using OrderedCollections
import JSON

# import Python matplotlib
plt = pyimport("matplotlib.pyplot")


"""
An internal utility fn to count only the number of defects (excludes start and end nodes)
"""
function _count_defects(all_pts::OrderedDict)
    return length([ki for ki in keys(all_pts) if !(ki in ("start", "end"))])
end


function capture_start(all_pts::OrderedDict, ax)
    pt = plt.ginput(1)[1]
    all_pts["start"] = [pt[1], pt[2]]

    ax.plot(pt[1], pt[2], "r^")
    ax.text(pt[1], pt[2], "start")

    return all_pts
end


function capture_end(all_pts::OrderedDict, ax)
    pt = plt.ginput(1)[1]
    all_pts["end"] = [pt[1], pt[2]]

    ax.plot(pt[1], pt[2], "rv")
    ax.text(pt[1], pt[2], "end")

    return all_pts
end


function capture_defect(all_pts::OrderedDict, ax; how_many=1)
    for _ in 1:how_many
        pt = plt.ginput(1)[1]
        n_defects = _count_defects(all_pts)
        all_pts["$(n_defects + 1)"] = [pt[1], pt[2]]

        ax.plot(pt[1], pt[2], "k.")
        ax.text(pt[1], pt[2], "$(n_defects + 1)")
    end

    return all_pts
end


function capture_segment(all_pts::OrderedDict, compulsory_edges::Array, ax)
    xs = []
    ys = []
    for _ in 1:2
        pt = plt.ginput(1)[1]
        
        n_defects = _count_defects(all_pts)
        all_pts["$(n_defects + 1)"] = [pt[1], pt[2]]

        push!(xs, pt[1])
        push!(ys, pt[2])

        ax.plot(pt[1], pt[2], "r.")
        ax.text(pt[1], pt[2], "$(n_defects + 1)")
    end

    # append edge node names for compulsory edges
    n_defects = _count_defects(all_pts)
    push!(compulsory_edges, ["$(n_defects - 1)", "$n_defects"])

    # plot connecting line + direction
    v1 = [xs[1], ys[1]]
    v2 = [xs[2], ys[2]]
    dv = v2 - v1


    ax.plot(xs, ys, "r-")
    ax.annotate("", xy=(v1 + 0.2*dv), xytext=v1, arrowprops=Dict("arrowstyle"=>"->", "color"=>"k"))

    return all_pts, compulsory_edges
end


"""
Creates a blank matplotlib canvas for use with ginput
"""
function make_blank_canvas()
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 10])
    ax.grid(linewidth=0.2)
    ax.set_title("Linear task network designer")

    plt.show()
    return fig, ax
end


function main_loop(all_pts, compulsory_edges, ax)
    start_node_specified = false
    end_node_specified = false
    
    # enter state loop
    while true
        print("\n--> (s)tart, (e)nd, (c)compulsory, (d)efect, (q)uit, sa(v)e: ")
        user_input = readline()

        if user_input == "s"
            println("start node selected")
            if start_node_specified
                println("a start node has already been specified")
                continue
            end
            all_pts = capture_start(all_pts, ax)
            start_node_specified = true
        elseif user_input == "e"
            println("end node selected")
            if end_node_specified
                println("an end node has already been specified")
                continue
            end
            all_pts = capture_end(all_pts, ax)
            end_node_specified = true
        elseif user_input == "c"
            println("compulsory segment selected")
            all_pts, compulsory_edges = capture_segment(all_pts, compulsory_edges, ax)
        elseif user_input == "d"
            print("defect node selected. How many defects? ")
            how_many = parse(Int, readline())
            all_pts = capture_defect(all_pts, ax, how_many=how_many)
        elseif user_input == "q"
            println("exit selected")
            break
        elseif user_input == "v"
            println("saving as json")
            print("provide a filename (e.g. test_network.json): ")
            save_filename = readline()

            # write to disk
            all_data = Dict(
                "points" => all_pts,
                "compulsory_edges" => compulsory_edges
            )

            open(save_filename, "w") do f
                JSON.print(f, all_data)
            end
        else
            println("unknown input, try again")
        end
    end

    return all_pts, compulsory_edges
end

# -------------------------------------------------------------------------------------------------

# data setup
# variables are defined here so that they can be accessed within the module's global scope
all_pts = OrderedDict{String, Vector{Float64}}()
compulsory_edges = []

# establish canvas
fig, ax = make_blank_canvas()

# main user input loop
all_pts, compulsory_edges = main_loop(all_pts, compulsory_edges, ax)

end