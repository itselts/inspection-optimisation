using CSV
using DataFrames
using JuMP
using Pipe
using Gurobi
using PyCall
using Random
using OrderedCollections
using LinearAlgebra
using ArgParse
using Parquet 
using Dates
# Scenario 1: Using the REPL model 
# Run optimise_ws.jl to populate the REPL with the relevant variables. To run optimise_ws.jl, uncomment all the arg parse stuff.
# That will solve the warm start and populate the variables, while also ingesting the GetOptimiseRunInspectionInputItems (Database values)
# The warm start solve will be scaled by 100

request_ID = "6D024294-5BF8-444C-86EC-6817F13766F4"
s_date_str = "2023-06-23"
s_time_str = "07:00"
DATE = DateTime(s_date_str)
DOW = dayname(DATE)
TIME = Time(s_time_str)


##### PLOTTING SOLUTION FUNCTION #####
plt = PyCall.pyimport("matplotlib.pyplot")

function plot_soln(x)
    fig = plt.figure()
    ax = fig.gca()
    ax.grid(linewidth=0.2)
    ax.set_title("Tour solution")

    # Plotting all inspections
    ax.plot(vcat(parse.(Float64, seg_df[:, "LongStart"]), parse.(Float64, seg_df[:, "LongEnd"])), vcat(parse.(Float64, seg_df[:, "LatStart"]), parse.(Float64, seg_df[:, "LatEnd"])), "k.")

    # Plotting all defects
    for row in eachrow(def_df)
        ax.plot(parse(Float64, row["LongStart"]), parse(Float64, row["LatStart"]), "yv", markersize = 5)
    end

    # Plot depot_df
    for row in eachrow(depot_df)
        ax.plot(parse(Float64, row["Long"]), parse(Float64, row["Lat"]), "m^", markersize=15)
    end

    # Plotting solution
    for k in CREWS
        for p1 in NODES, p2 in NODES
            if round(Int, value(x[p1, p2, k])) == 1
                sx, sy = NODES_DATA[p1]
                ex, ey = NODES_DATA[p2]
                
                # Add vectors
                dp = [ex-sx, ey-sy]
                ax.annotate("", xy=([sx, sy] + 0.2*dp), xytext=[sx ,sy], arrowprops=Dict("arrowstyle"=>"->", "color"=>"g"))
                # ax.annotate(string(value(t[p2]))[1:5], xy=[ex, ey], xytext=[ex ,ey])

                if (length(collect(eachsplit(p1, "_"))) == 4) & (length(collect(eachsplit(p2, "_"))) == 4)
                    if collect(eachsplit(p1, "_"))[1:3] == collect(eachsplit(p2, "_"))[1:3]
                        ax.plot([sx, ex], [sy, ey], linestyle="-", color="r", linewidth=1)
                    else
                        ax.plot([sx, ex], [sy, ey], linestyle="-", color="b", linewidth=1)
                    end
                elseif (length(collect(eachsplit(p1, "_"))) == 1) || (length(collect(eachsplit(p2, "_"))) == 1)
                    ax.plot([sx, ex], [sy, ey], linestyle="-", color="y", linewidth=1)
                else
                    ax.plot([sx, ex], [sy, ey], linestyle="-", color="b", linewidth=1)
                end
            end
        end
    end
    if isdir(string("../runs/", request_ID, "/results/solutions_plot.png"))
        rm(string("../runs/", request_ID, "/results/solutions_plot.png"), recursive=false)
    end
    plt.savefig(string("../runs/", request_ID, "/results/solutions_plot.png"))
    #plt.show()

    println("Plots saved.")
end

###### DATA INGESTION #####
# inspection segments
seg_df = @pipe read_parquet(string("../runs/", request_ID , "/data/v_GetOptimiseInspectionInputItems.parquet"))  |> DataFrame

# defects
def_df = @pipe read_parquet(string("../runs/", request_ID , "/data/v_GetOptimiseInputItems.parquet"))  |> DataFrame

# depots
depot_df = @pipe read_parquet(string("../runs/", request_ID , "/data/v_DomainValueDepot.parquet"))  |> DataFrame

# crews
crew_df = @pipe read_parquet(string("../runs/", request_ID , "/data/v_GetOptimiseInputCrews.parquet"))  |> DataFrame

# dist/time matrix
pairwise_cache = @pipe read_parquet(string("../runs/", request_ID , "/data/v_GetDistTimeCacheForOptRun.parquet"))  |> DataFrame # cache pairwise
pairwise_cache = pairwise_cache[completecases(pairwise_cache), :]

# pairwise_cache[[!x for x in completecases(pairwise_cache)], [:PointAID_FK, :PointBID_FK, :Distance, :Time]] # Find all rows with missing entries

pairwise1 = copy(pairwise_cache[!, [:PointAID_FK, :PointBID_FK, :Distance, :Time]])
pairwise1.PointAID_FK = string.(pairwise1.PointAID_FK) # Type conversions
pairwise1.PointBID_FK = string.(pairwise1.PointBID_FK)
pairwise1.Time = Float64.(pairwise1.Time)
pairwise1.Distance = Float64.(pairwise1.Distance)

pairwise2 = seg_df[!, [:SegStartIdentifier, :SegEndIdentifier, :EstQty, :EstDuration]] # segments pairwise
pairwise2.EstDuration = Float64.(pairwise2.EstDuration)  # Type conversions
pairwise2.EstQty = Float64.(pairwise2.EstQty) ./ 1000 # Metres to kilometres
rename!(pairwise2, :SegStartIdentifier => :PointAID_FK, :SegEndIdentifier => :PointBID_FK, :EstQty => :Distance, :EstDuration => :Time)

pairwise = vcat(pairwise1, pairwise2) # Full pairwise list 

dist_df = @pipe unstack(pairwise, :PointAID_FK, :PointBID_FK, :Distance) # Pairwise to DataFrame 
time_df = @pipe unstack(pairwise, :PointAID_FK, :PointBID_FK, :Time)


ORIGINS = string.(dist_df[:, :PointAID_FK])
DESTINATIONS = names(select(dist_df, Not(:PointAID_FK)))

dist_df = select(dist_df, Not(:PointAID_FK)) # Remove the row index column from unstack (Index is stored in ORIGINS)
time_df = select(time_df, Not(:PointAID_FK))

dist_mat = @pipe dist_df |> Matrix |> JuMP.Containers.DenseAxisArray(_, ORIGINS, DESTINATIONS) # DataFrame to DenseAxisArray
time_mat = @pipe time_df |> Matrix |> JuMP.Containers.DenseAxisArray(_, ORIGINS, DESTINATIONS)


# CREW START AND END 
crew_start_end = Dict{}()
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    for i in 1:crew_num
        depot_code_start = crew_row["StartDepotCode"]
        depot_code_end = crew_row["EndDepotCode"]
        
        depot_refnum_start = depot_df[depot_df.dv_code .== depot_code_start, "DepotRefNumber"][1]
        depot_refnum_end = depot_df[depot_df.dv_code .== depot_code_end, "DepotRefNumber"][1]
        
        crew_start_end[string(crew_row["CrewTypeCode"], "_", i)] = Dict("start" => string(depot_refnum_start), "end" => string(depot_refnum_end))
    end
end

# CREW SHIFT TIMES
SHIFT_END = Dict{String, Int64}()
SHIFT_START = Dict{String, Int64}()

for row in eachrow(crew_df)
    crew_num = row[string(DOW, "NumCrew")]
    for i in 1:crew_num
        crew = string(row["CrewTypeCode"], "_", i)
        start_time = Time(row[string(DOW, "StartTime")] - Dates.Hour(12))
        end_time = Time(row[string(DOW, "EndTime")] - Dates.Hour(12))

        SHIFT_START[crew] = 60*Hour(start_time).value + Minute(start_time).value
        SHIFT_END[crew] = 60*Hour(end_time).value + Minute(end_time).value
    end
end


##### SETS #####
NODES_DATA = Dict()
for row in eachrow(seg_df)
    NODES_DATA[string(row["SegStartIdentifier"])] = (parse(Float64, row["LongStart"]), parse(Float64, row["LatStart"]))
    NODES_DATA[string(row["SegEndIdentifier"])] = (parse(Float64, row["LongEnd"]), parse(Float64, row["LatEnd"]))
end

for row in eachrow(def_df)
    NODES_DATA[string(row["ItemIdentifier"], "0")] =  (parse(Float64, row["LongStart"]), parse(Float64, row["LatStart"])) # The START of each defect
end

for row in eachrow(depot_df)
    NODES_DATA[string(row["DepotRefNumber"])] = (parse(Float64, row["Long"]), parse(Float64, row["Lat"]))
end

NODES = string.(collect(keys(NODES_DATA)))
SEGMENTS = [string(x[1]) for x in eachrow(seg_df[!, "ItemIdentifier"] )]
SEGMENT_POINTS = []

for row in eachrow(seg_df)
    push!(SEGMENT_POINTS, row["SegStartIdentifier"])
    push!(SEGMENT_POINTS, row["SegEndIdentifier"])
end

DEFECTS = [join([row["ItemIdentifier"], "0"]) for row in eachrow(def_df)] # Defect START points
JOBS = string.(vcat(SEGMENT_POINTS, DEFECTS))
DEPOTS = string.(depot_df[!, "DepotRefNumber"])
CREW_TYPES = string.(Set(crew_df[!, :CrewTypeCode]))

CREWS = []
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    for i in 1:crew_num
        push!(CREWS, string(crew_row["CrewTypeCode"], "_", i))
    end
end

CREW_CAPABILITY = Dict(crew => Array{String}(undef,0) for crew in CREWS) # List of
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    for i in 1:crew_num
        crew_ID = string(crew_row["CrewTypeCode"], "_", i)
        crew_type = crew_row["CrewTypeCode"]
        for seg_row in eachrow(seg_df)
            if occursin(crew_type, seg_row["CrewType_FKs"])
                push!(CREW_CAPABILITY[crew_ID], string(seg_row["ItemIdentifier"], "0")) # Start of segment
                push!(CREW_CAPABILITY[crew_ID], string(seg_row["ItemIdentifier"], "1")) # End of segment
            end
        end

        for def_row in eachrow(def_df)
            if occursin(crew_type, def_row["CrewType_FKs"])
                push!(CREW_CAPABILITY[crew_ID], string(def_row["ItemIdentifier"], "0")) # Start of segment
                push!(CREW_CAPABILITY[crew_ID], string(def_row["ItemIdentifier"], "1")) # End of segment
            end
        end
    end
end


##### CONSTANTS/MODEL INPUTS #####
REVERSIBLE = Dict(string(row["ItemIdentifier"]) => row["Reversible"] for row in eachrow(seg_df)) # Dictionary comprehension

DUR = Dict{String, Int}() # Estimated DURation of each defect
for row in eachrow(def_df)
    DUR[string(row["ItemIdentifier"], "0")] = row["EstDuration"]
end

SCORES = @pipe CSV.File(string("../runs/", request_ID , "/outputs/scores.csv"))  |> DataFrame # Job scores
SCORES = Dict{String, Float64}(string(row["ItemIdentifier"], "0") => row["Score"] for row in eachrow(SCORES)) # Append 0 to ID for start of defect


##### MODELLING #####
model = Model(Gurobi.Optimizer)

@variable(model, x[NODES, NODES, CREWS], Bin)
@variable(model, y[JOBS, CREWS], Bin)
@variable(model, t[JOBS] >= 0)

M = 2000

# 0. Graph reduction from preprocessing logic
for (k, v) in zip(keys(time_mat), values(time_mat))
    if ismissing(v)
        local i, j = k[1], k[2]
        for c in CREWS
            try 
                @constraint(model, x[i, j, c] == 0)
            catch ex
                continue
            end
        end
    end
end

# 1. inspection segments must be done, as well as if it is bidirectional
for s in SEGMENTS
    local i = join([s, "0"]) # Start point of segment is segment ID concatenated with 0 
    local j = join([s, "1"]) # End point of segment is segment ID concatenated with 1
    
    if REVERSIBLE[s] == 1
        @constraint(model, sum(x[i, j, k] for k in CREWS) + sum(x[j, i, k] for k in CREWS) == 1)
    else
        @constraint(model, sum(x[i, j, k] for k in CREWS) == 1)
        @constraint(model, sum(x[j, i, k] for k in CREWS) == 0)
        #@constraint(model, sum(x[i, j_other, k] for k in CREWS, j_other in setdiff(NODES, [j])) == 0)
    end
    
end

# 2. conservation of flow. If a crew enters, crew must also leave 
for j in NODES
    for k in CREWS
        @constraint(model, sum(x[i, j, k] for i in NODES) == sum(x[j, i, k] for i in NODES))
    end
end

# 2. at most one crew can be active leaving jobs
for i in JOBS
    @constraint(model, sum(x[i, j, k] for j in NODES, k in CREWS) <= 1)
end

# 3. at most one crew can be active entering jobs
for j in JOBS
    @constraint(model, sum(x[i, j, k] for i in NODES, k in CREWS) <= 1)
end

# 4. If job j is active, one of the entering edges attached to j must be active (linking)
for j in JOBS, k in CREWS
    @constraint(model, y[j, k] == sum(x[i, j, k] for i in NODES))
end

# 5. If job i is active, one of the exiting edges must be active (linking)
for i in JOBS, k in CREWS
    @constraint(model, y[i, k] == sum(x[i, j, k] for j in NODES))
end

# 6. each vehicle leaving and returning
for k in CREWS
    local START = crew_start_end[k]["start"]
    local END = crew_start_end[k]["end"]
    @constraint(model, sum(x[START, j, k] for j in NODES) == 1)
    @constraint(model, sum(x[i, END, k] for i in NODES) == 1)
end

# 7. Can't enter or leave a depot that is not the crews designated depot
for k in CREWS
    local START = crew_start_end[k]["start"]
    local END = crew_start_end[k]["end"]
    for OTHER_DEPOT in setdiff(setdiff(DEPOTS, [START]), [END])
        @constraint(model, sum(x[i, OTHER_DEPOT, k] for i in NODES) == 0) # Can't enter another depot
        @constraint(model, sum(x[OTHER_DEPOT, j, k] for j in NODES) == 0) # Can't leave another depot
    end
end

# 8. cannot travel to self or depot
for k in CREWS
    local START = crew_start_end[k]["start"]
    local END = crew_start_end[k]["end"]
    @constraint(model, sum(x[i, i, k] for i in NODES) == 0)
    @constraint(model, x[START, END, k] == 0)
    @constraint(model, x[END, START, k] == 0)
end

# 9. time constraints (subtours)
for i in JOBS, j in JOBS
    if j in DEFECTS
        try
            @constraint(model, t[j] >= t[i] + time_mat[i,j] + DUR[j] - M*(1 - sum(x[i,j,k] for k in CREWS)))            
        catch ex
        end
    else
        try
            @constraint(model, t[j] >= t[i] + time_mat[i,j] - M*(1 - sum(x[i,j,k] for k in CREWS)))
        catch ex
        end
    end
    
end

# 10. shift start and end time
for j in JOBS
    if j in DEFECTS
        try # If time matrix entry is missing
            @constraint(model, t[j] >= sum(x[crew_start_end[k]["start"], j, k]*(SHIFT_START[k] + DUR[j] + time_mat[crew_start_end[k]["start"], j]) for k in CREWS)) # Start time
        catch ex
        end
        try
            begin @constraint(model, t[j] <= sum(x[j, crew_start_end[k]["end"], k]*(SHIFT_END[k] - time_mat[crew_start_end[k]["end"], j]) for k in CREWS) 
            + M*(1 - sum(x[j, crew_start_end[k]["end"], k] for k in CREWS))) end # End time
        catch ex
        end
    else    
        try # If time matrix entry is missing
            @constraint(model, t[j] >= sum(x[crew_start_end[k]["start"], j, k]*(SHIFT_START[k] + time_mat[crew_start_end[k]["start"], j]) for k in CREWS)) # Start time
        catch ex
        end
        try
            begin @constraint(model, t[j] <= sum(x[j, crew_start_end[k]["end"], k]*(SHIFT_END[k] - time_mat[crew_start_end[k]["end"], j]) for k in CREWS) 
            + M*(1 - sum(x[j, crew_start_end[k]["end"], k] for k in CREWS))) end # End time
        catch ex
        end
    end
end

# 11. Crew type constraints
for i in JOBS
    for k in CREWS
        if !(i in CREW_CAPABILITY[k])
            @constraint(model, sum(x[i, j, k] for j in NODES) == 0)
            @constraint(model, sum(x[j, i, k] for j in NODES) == 0)
        end
    end
end


# objective
@expression(model, tot_time, sum(coalesce(time_mat[i, j], 10000) * x[i, j, k] for i in NODES, j in NODES, k in CREWS))
@expression(model, tot_score, sum(SCORES[i]*y[i, k] for i in DEFECTS, k in CREWS))
@expression(model, tot_defects, sum(y[i, k] for i in DEFECTS, k in CREWS))
@expression(model, tot_dist, sum(coalesce(dist_mat[i, j], 10000) * x[i, j, k] for i in NODES, j in NODES, k in CREWS))
@objective(model, Max, tot_score - tot_time)






##### WARM START FULL SOLUTION FROM PRIOR SOLVE #####
warm_start_ID = "DB00491A-36EB-49E6-9295-04C7E52A6AD6"
warm_start = CSV.read(string("../runs/", warm_start_ID , "/warm_start/variables/all_variables.csv"), DataFrame)


#=
##### INGESTING FULL WARM START #####
for row in eachrow(warm_start)
    local var, val = row["Variable"], row["Value"]
    if occursin("x[", var)
        local i, j, k = @pipe var |> split(_, "[")[2] |> split(_, "]")[1] |> split(_, ",")
        set_start_value(x[i,j,k], val)
    elseif occursin("y[", var)
        i, k = @pipe var |> split(_, "[")[2] |> split(_, "]")[1] |> split(_, ",")
        set_start_value(y[i, k], val)
    else
        i = @pipe var |> split(_, "[")[2] |> split(_, "]")[1]
        set_start_value(t[i], val) 
    end
end
=#


##### INGESTING AND PROCESSING THE WARM START FROM THE DATABASE #####
set_start_value.(all_variables(model), 0)

seg_df.DefaultShiftStartOffset = Float64.(seg_df.DefaultShiftStartOffset) # Type conversion from decimal to float
seg_df.DefaultSequence = Float64.(seg_df.DefaultSequence)

for row in eachrow(crew_df) 
    crew_num = row[string(DOW, "NumCrew")]
    for i in 1:crew_num
        crew = string(row["CrewTypeCode"], "_", i)
        crew_type = row["CrewTypeCode"]
        local start = crew_start_end[crew]["start"]
        seg_count = maximum(seg_df[seg_df.CrewType_FKs .== crew_type, :DefaultSequence]) # Number of segments for 1 inspection route

        mask = (seg_df.DefaultSequence .== 1) .&& (seg_df.CrewType_FKs .== crew_type)
        local to = string(seg_df[mask, "ItemIdentifier"][1]) # To segment ID

        set_start_value(x[start, string(to, "0"), crew], 1) # Start depot to first segment start point
        set_start_value(t[string(to, "0")], SHIFT_START[crew] + seg_df[mask, "DefaultShiftStartOffset"][1]) # First segment start time
        set_start_value(t[string(to, "1")], SHIFT_START[crew] + seg_df[mask, "DefaultShiftStartOffset"][1] + seg_df[mask, "EstDuration"][1]) # First segment end time

        local j = 2
        while j <= seg_count
            mask = (seg_df.DefaultSequence .== j) .&& (seg_df.CrewType_FKs .== crew_type)
            local from = to # From segment ID
            to = string(seg_df[mask, "ItemIdentifier"][1]) # To segment ID

            set_start_value(x[string(from, "0"), string(from, "1"), crew], 1) # Inspection segment
            set_start_value(x[string(from, "1"), string(to, "0"), crew], 1) # From segment end, to another segment start

            set_start_value(y[string(from, "0"), crew], 1) # Job start indicator
            set_start_value(y[string(from, "1"), crew], 1) # Job finish indicator

            set_start_value(t[string(to, "0")], SHIFT_START[crew] + seg_df[mask, "DefaultShiftStartOffset"][1]) # Segment start time 
            set_start_value(t[string(to, "1")], SHIFT_START[crew] + seg_df[mask, "DefaultShiftStartOffset"][1] + seg_df[mask, "EstDuration"][1]) # Segment finish time 
            
            j += 1
        end
        set_start_value(x[string(to, "0"), string(to, "1"), crew], 1) # Last segment

        set_start_value(y[string(to, "0"), crew], 1) # Last segment indicator
        set_start_value(y[string(to, "1"), crew], 1) # Last segment indicator

        finish = crew_start_end[crew]["end"]
        set_start_value(x[string(to, "1"), finish, crew], 1) # Last segment end to depot finish
    end
end

X = all_variables(model)
processed_warm_start = DataFrame(zip(string.(X), start_value.(X))) # Warm start in database processed
rename!(processed_warm_start, ["Variable", "Value"])



##### COMPARE #####
for row in eachrow(warm_start)
    var = row["Variable"]
    value = row["Value"]
    if !(processed_warm_start[processed_warm_start.Variable .== var, "Value"][1] == value)
        if occursin("t[", var)
            if abs(value - processed_warm_start[processed_warm_start.Variable .== var, "Value"][1]) >= 0.001
                print(var, ": ")
                println("Warm start = ", value, " Database = ", processed_warm_start[processed_warm_start.Variable .== var, "Value"][1])
            end
        else
        print(var, ": ")
        println("Warm start = ", value, " Database = ", processed_warm_start[processed_warm_start.Variable .== var, "Value"][1])
        end
    end
end


##### SOLVING FOR THE INGESTED WARM START #####
set_optimizer_attribute(model, "SolutionLimit", 1)
set_optimizer_attribute(model, "MIPFocus", 1)
optimize!(model)

println("***************************************")

##### SOLVING WITH FULL WARM START #####
X = all_variables(model)
x_solution = value.(X)
set_start_value.(X, x_solution)

MOIU.reset_optimizer(model) # Resets attributes?
set_optimizer_attribute(model, "SolutionLimit", 5)
set_optimizer_attribute(model, "TimeLimit", 30)
set_optimizer_attribute(model, "NoRelHeurTime", 30)
optimize!(model)

