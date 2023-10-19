# Activating the project environment
using Pkg
println("----------------------------------------")
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
    Pkg.instantiate()
end

using CSV
using DataFrames
using JuMP
using Pipe
using Gurobi
using Random
using OrderedCollections
using LinearAlgebra
using ArgParse
using Parquet 
using Dates

parser = ArgParseSettings()
@add_arg_table parser begin
    "--request_ID"
        arg_type = String
        required = true
    "--s_date"
        arg_type = String
        required = true
    "--s_time"
        arg_type = String
        required = true
end

args = parse_args(parser)
request_ID = string(args["request_ID"])
s_date_str = string(args["s_date"])
s_time_str = string(args["s_time"])
DATE = DateTime(s_date_str)
DOW = dayname(DATE)
TIME = Time(s_time_str)

# request_ID = "EF4DD65B-B09B-47DB-ABF4-2EB8D9B2A211"
# s_date_str = "2023-06-23"
# s_time_str = "07:00"
# DATE = DateTime(s_date_str)
# DOW = dayname(DATE)
# TIME = Time(s_time_str)

###### DATA INGESTION #####
print("Ingesting and transforming data...")

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

"Function that checks if any pair has a missing time or distance value."
function pairwise_check(pairwise_cache)
    if size(pairwise_cache) != size(pairwise_cache[completecases(pairwise_cache), :])
        println("")
        println("******** ROWS IN THE CACHE IS MISSING TIME OR DIST VALUE ********")
        println(size(pairwise_cache)[1] - size(pairwise_cache[completecases(pairwise_cache), :])[1], " rows missing.")
        println("")
    end
end
pairwise_check(pairwise_cache)

pairwise_cache = pairwise_cache[completecases(pairwise_cache), :]

pairwise1 = copy(pairwise_cache[!, [:PointAID_FK, :PointBID_FK, :Distance, :Time]]) # defects pairwise
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

# NODES FROM THE CACHE TABLE (c.f TO NODES FROM THE INPUT ITEM TABLES)
ORIGINS = string.(dist_df[:, :PointAID_FK]) # When unstacking, first column is PointAID_FK, which are the origins.
DESTINATIONS = names(select(dist_df, Not(:PointAID_FK))) # Get the column names excluding POINTAID_FK, which are the destinations.

dist_df = select(dist_df, Not(:PointAID_FK)) # Remove the row index column from unstack (Index is stored in ORIGINS)
time_df = select(time_df, Not(:PointAID_FK))

dist_mat = @pipe dist_df |> Matrix |> JuMP.Containers.DenseAxisArray(_, ORIGINS, DESTINATIONS) # DataFrame to DenseAxisArray
time_mat = @pipe time_df |> Matrix |> JuMP.Containers.DenseAxisArray(_, ORIGINS, DESTINATIONS)


# CREW START AND END 
CREW_START_END = Dict{}()
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    for i in 1:crew_num
        depot_code_start = crew_row["StartDepotCode"]
        depot_code_end = crew_row["EndDepotCode"]
        
        depot_refnum_start = depot_df[depot_df.dv_code .== depot_code_start, "DepotRefNumber"][1]
        depot_refnum_end = depot_df[depot_df.dv_code .== depot_code_end, "DepotRefNumber"][1]
        
        CREW_START_END[string(crew_row["CrewTypeCode"], "_", i)] = Dict("start" => string(depot_refnum_start), "end" => string(depot_refnum_end))
    end
end


# CREW SHIFT TIMES
SHIFT_END = Dict{String, Int64}()
SHIFT_START = Dict{String, Int64}()

for row in eachrow(crew_df)
    crew_num = row[string(DOW, "NumCrew")]
    if crew_num == 0
        println("********** There are no crews rostered on ", DOW, " for ", row["CrewTypeCode"] ,". **********")
    end
    for i in 1:crew_num
        crew = string(row["CrewTypeCode"], "_", i)
        start_time = Time(row[string(DOW, "StartTime")] - Dates.Hour(12))
        end_time = Time(row[string(DOW, "EndTime")] - Dates.Hour(12))

        SHIFT_START[crew] = 60*Hour(start_time).value + Minute(start_time).value
        SHIFT_END[crew] = 60*Hour(end_time).value + Minute(end_time).value
    end
end


##### SETS FROM INPUT DATA #####
NODES_DATA = Dict() # NODES FROM THE INPUT ITEM TABLES (AS OPPOSED TO NODES FROM THE CACHE TABLE)
for row in eachrow(seg_df)
    NODES_DATA[string(row["SegStartIdentifier"])] = (parse(Float64, row["LongStart"]), parse(Float64, row["LatStart"]))
    NODES_DATA[string(row["SegEndIdentifier"])] = (parse(Float64, row["LongEnd"]), parse(Float64, row["LatEnd"]))
end

for row in eachrow(def_df)
    NODES_DATA[string(row["ItemIdentifier"], "0")] = (parse(Float64, row["LongStart"]), parse(Float64, row["LatStart"])) # The START of each defect
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

"Function that checks for points that are in the input items, but not in the cache."
function cache_getoptitems_check()
    cache = setdiff(DESTINATIONS, DEPOTS)
    if Set(cache) != Set(JOBS)
        println("******** THERE ARE POINTS IN EITHER v_GetOptimiseInspectionInputItems OR v_GetOptimiseInputItems BUT ARE NOT IN v_GetDistTimeCacheForOptRun ******")
        for point in JOBS
            if !(point in setdiff(DESTINATIONS, DEPOTS))
                println("Point ", point, " is missing from the cache.")
            end
        end
    end
end
cache_getoptitems_check()

##### CONSTANTS/MODEL INPUTS #####
REVERSIBLE = Dict(string(row["ItemIdentifier"]) => row["Reversible"] for row in eachrow(seg_df)) # Dictionary comprehension

DUR = Dict{String, Int}() # Estimated DURation of each defect
for row in eachrow(def_df)
    DUR[string(row["ItemIdentifier"], "0")] = row["EstDuration"]
end

SCORES = @pipe CSV.File(string("../runs/", request_ID , "/outputs/scores.csv"))  |> DataFrame # Job scores
SCORES = Dict{String, Float64}(string(row["ItemIdentifier"], "0") => row["Score"] for row in eachrow(SCORES)) # Append 0 to ID for start of defect

println("ok")

##### MODELLING #####
model = Model(Gurobi.Optimizer)

print("Creating the model...")

@variable(model, x[NODES, NODES, CREWS], Bin)
@variable(model, y[JOBS, CREWS], Bin)
@variable(model, t[JOBS] >= 0)

# Big M constants
M1_1 = maximum(values(SHIFT_END)) + maximum(values(DUR)) + maximum(skipmissing(time_mat.data))
M1_2 = maximum(values(SHIFT_END)) + maximum(skipmissing(time_mat.data))
M2 = maximum(values(SHIFT_END)) 

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
    local START = CREW_START_END[k]["start"]
    local END = CREW_START_END[k]["end"]
    @constraint(model, sum(x[START, j, k] for j in NODES) <= 1)
    @constraint(model, sum(x[i, END, k] for i in NODES) <= 1)
end

# 7. Can't enter or leave a depot that is not the crews designated depot
for k in CREWS
    local START = CREW_START_END[k]["start"]
    local END = CREW_START_END[k]["end"]
    for OTHER_DEPOT in setdiff(setdiff(DEPOTS, [START]), [END])
        @constraint(model, sum(x[i, OTHER_DEPOT, k] for i in NODES) == 0) # Can't enter another depot
        @constraint(model, sum(x[OTHER_DEPOT, j, k] for j in NODES) == 0) # Can't leave another depot
    end
end

# 8. cannot travel to self or depot
for k in CREWS
    local START = CREW_START_END[k]["start"]
    local END = CREW_START_END[k]["end"]
    @constraint(model, sum(x[i, i, k] for i in NODES) == 0)
    @constraint(model, x[START, END, k] == 0)
    @constraint(model, x[END, START, k] == 0)
end

# 9. time constraints (subtours)
for i in JOBS, j in JOBS
    if j in DEFECTS
        try
            @constraint(model, t[j] >= t[i] + time_mat[i,j] + DUR[j] - M1_1*(1 - sum(x[i,j,k] for k in CREWS)))            
        catch ex
        end
    else
        try
            @constraint(model, t[j] >= t[i] + time_mat[i,j] - M1_2*(1 - sum(x[i,j,k] for k in CREWS))) 
        catch ex
        end
    end
end

# 10. shift start and end time
for j in JOBS
    if j in DEFECTS
        try # If time matrix entry is missing
            @constraint(model, t[j] >= sum(x[CREW_START_END[k]["start"], j, k]*(SHIFT_START[k] + DUR[j] + time_mat[CREW_START_END[k]["start"], j]) for k in CREWS)) # Start time
        catch ex
        end
        try
            begin @constraint(model, t[j] <= sum(x[j, CREW_START_END[k]["end"], k]*(SHIFT_END[k] - time_mat[CREW_START_END[k]["end"], j]) for k in CREWS) 
            + M2*(1 - sum(x[j, CREW_START_END[k]["end"], k] for k in CREWS))) end # End time
        catch ex
        end
    else    
        try # If time matrix entry is missing
            @constraint(model, t[j] >= sum(x[CREW_START_END[k]["start"], j, k]*(SHIFT_START[k] + time_mat[CREW_START_END[k]["start"], j]) for k in CREWS)) # Start time
        catch ex
        end
        try
            begin @constraint(model, t[j] <= sum(x[j, CREW_START_END[k]["end"], k]*(SHIFT_END[k] - time_mat[CREW_START_END[k]["end"], j]) for k in CREWS) 
            + M2*(1 - sum(x[j, CREW_START_END[k]["end"], k] for k in CREWS))) end # End time
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

println("ok")

##### WARM START STORED IN DATABASE #####
if !isempty(seg_df) # Check if it is a pure defect optimisation with no inspections
    print("Processing the warm starts for the inspection runs...")
    set_start_value.(all_variables(model), 0)

    seg_df.DefaultShiftStartOffset = Float64.(seg_df.DefaultShiftStartOffset) # Type conversion from decimal to float
    seg_df.DefaultSequence = Float64.(seg_df.DefaultSequence)
    
    for row in eachrow(crew_df) 
        crew_num = row[string(DOW, "NumCrew")]
        for i in 1:crew_num
            crew = string(row["CrewTypeCode"], "_", i)
            crew_type = row["CrewTypeCode"]
            local start = CREW_START_END[crew]["start"]
            try # Check to see if an inspection crew. If it is, ingest warm start
                seg_count = maximum(seg_df[seg_df.CrewType_FKs .== crew_type, :DefaultSequence]) # Number of segments for 1 inspection route
                mask = (seg_df.DefaultSequence .== 1) .&& (seg_df.CrewType_FKs .== crew_type) # Mask to the first inspection
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
        
                finish = CREW_START_END[crew]["end"]
                set_start_value(x[string(to, "1"), finish, crew], 1) # Last segment end to depot finish

            catch ex
                continue
            end
        end
    end
    println("ok")

    ##### SOLVING FOR THE INGESTED WARM START #####
    println("Readjusting the warm start...\n")
    set_optimizer_attribute(model, "SolutionLimit", 1)
    set_optimizer_attribute(model, "MIPFocus", 1)
    optimize!(model)

    ##### SETTING START VARIABLES WITH FULL WARM START #####
    X = all_variables(model)
    x_solution = value.(X)
    set_start_value.(X, x_solution)

    println("\n\n\nWarm start successfully processed.")
end

#=
##### WARM START FROM PRIOR SOLVE #####
set_start_value.(all_variables(model), 0)
warm_start = CSV.read(string("../runs/", request_ID , "/warm_start/variables/all_variables.csv"), DataFrame)

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
        set_start_value(t[i], val/100) # Warm start scales the time by 100
    end
end
=#

# write mps model file 
#write_to_file(model, "../example2/results/model.mps")

##### NO REL SOLVE #####
println("Beginning phase 2 NoRel heuristic solve..\n\n")
MOIU.reset_optimizer(model) # Resets attributes? Seems like this is broken
set_optimizer_attribute(model, "SolutionLimit", 10000)
set_optimizer_attribute(model, "TimeLimit", 1020)
set_optimizer_attribute(model, "NoRelHeurTime", 1020)
optimize!(model)


##### PROCESSING ORDER OF POINTS ######
SEQUENCE = Dict{String, Vector{String}}()
for k in CREWS
    i = 1
    start = CREW_START_END[k]["start"]
    if !(isnothing(findfirst([round(value(x[start, j, k])) for j in NODES] .!= 0)))
        SEQUENCE[k] = []
        to = @pipe findfirst([round(value(x[start, j, k])) for j in NODES] .!= 0) |> NODES[_] # Find the numeric index of the node that it leaves to, then get the node label.
        push!(SEQUENCE[k], start)
        push!(SEQUENCE[k], to)
        while to != CREW_START_END[k]["end"]
            from = to
            to = @pipe findfirst([round(value(x[from, j, k])) for j in NODES] .!= 0) |> NODES[_]
            push!(SEQUENCE[k], to)
        end
    end
end


# variable values
all_var_dict = Dict(zip(all_variables(model), value.(all_variables(model))))

##### OptRunMetricsOutput #####
metrics = Dict{String, Any}()
metrics["ModelIdentifier"] = "v1.0.0"
metrics["ResultStatusCode"] = primal_status(model)
metrics["SolutionFound"] = has_values(model)
metrics["TotalJobsCompleted"] = value(tot_defects)
metrics["TotalMissedJobs"] = length(NODES) - length(DEPOTS) -2*length(SEGMENTS) - value(tot_defects)
metrics["TotalPointsScore"] = objective_value(model)
metrics["TotalDistanceDriven"] = value(tot_dist)
metrics["RunTime"] = solve_time(model)
metrics["Gap"] = relative_gap(model)

println("OptRunOutputMetrics table created successfully.")

##### OptRunOutput #####
output = Dict{String, Any}()
itemtype = []
itemidentifier = []
optshift = []
crewtypeid = []
crewtypecode = []
crewnum = []
sequence = []
distancetravelforjob = []
pointscalculatedforjob = []
planstartdate = []
planenddate = []

# Defects
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    crew_type = crew_row["CrewTypeCode"]
    for i in 1:crew_num 
        crew = string(crew_row["CrewTypeCode"], "_", i)
        for def_row in eachrow(def_df)
            id = def_row["ItemIdentifier"]
            node = string(id, "0")
            if value(y[node, crew]) == 1
                push!(itemtype, def_row["ItemType"])
                push!(itemidentifier, def_row["ItemIdentifier"])
                push!(optshift, DATE+Minute(SHIFT_START[crew]))
                push!(crewtypeid, crew_row["CrewTypeID"])
                push!(crewtypecode, crew_row["CrewTypeCode"])
                push!(crewnum, i)
                seq_idx = findfirst(x -> x == node, SEQUENCE[crew])
                push!(sequence, seq_idx)
                push!(distancetravelforjob, dist_mat[SEQUENCE[crew][seq_idx-1], node])
                push!(pointscalculatedforjob, SCORES[node])
                push!(planstartdate, DATE + Minute(trunc(value(t[node]))) + Second(trunc(60*mod(value(t[node]), 1))) - Minute(trunc(DUR[node])) - Second(trunc(60*mod(DUR[node], 1))))
                push!(planenddate, DATE + Minute(trunc(value(t[node]))) + Second(trunc(60*mod(value(t[node]), 1))))
            end
        end
    end
end

# Segments
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    crew_type = crew_row["CrewTypeCode"]
    for i in 1:crew_num 
        crew = string(crew_row["CrewTypeCode"], "_", i)
        for seg_row in eachrow(seg_df)
            id = seg_row["ItemIdentifier"]
            start = string(seg_row["ItemIdentifier"], "0")
            finish = string(seg_row["ItemIdentifier"], "1")

            if value(y[start, crew]) == 1 # Segment starts
                push!(itemtype, "Segment Start")
                push!(itemidentifier, seg_row["ItemIdentifier"])
                push!(optshift, DATE+Minute(SHIFT_START[crew]))
                push!(crewtypeid, crew_row["CrewTypeID"])
                push!(crewtypecode, crew_row["CrewTypeCode"])
                push!(crewnum, i)
                seq_idx = findfirst(x -> x == start, SEQUENCE[crew])
                push!(sequence, seq_idx)
                push!(distancetravelforjob, dist_mat[SEQUENCE[crew][seq_idx-1], start])
                push!(pointscalculatedforjob, missing)
                push!(planstartdate, DATE + Minute(trunc(value(t[start]))) + Second(trunc(60*mod(value(t[start]), 1))))
                push!(planenddate, DATE + Minute(trunc(value(t[start]))) + Second(trunc(60*mod(value(t[start]), 1))))
            end

            if value(y[finish, crew]) == 1 # Segment ends
                push!(itemtype, "Segment End")
                push!(itemidentifier, seg_row["ItemIdentifier"])
                push!(optshift, DATE+Minute(SHIFT_START[crew]))
                push!(crewtypeid, crew_row["CrewTypeID"])
                push!(crewtypecode, crew_row["CrewTypeCode"])
                push!(crewnum, i)
                seq_idx = findfirst(x -> x == finish, SEQUENCE[crew])
                push!(sequence, seq_idx)
                push!(distancetravelforjob, dist_mat[SEQUENCE[crew][seq_idx-1], finish])
                push!(pointscalculatedforjob, missing)
                push!(planstartdate, DATE + Minute(trunc(value(t[finish]))) + Second(trunc(60*mod(value(t[finish]), 1))))
                push!(planenddate, DATE + Minute(trunc(value(t[finish]))) + + Second(trunc(60*mod(value(t[finish]), 1))))
            end
        end
    end
end

# Leave and return to depot
for crew_row in eachrow(crew_df)
    crew_num = crew_row[string(DOW, "NumCrew")]
    crew_type = crew_row["CrewTypeCode"]
    for i in 1:crew_num 
        crew = string(crew_row["CrewTypeCode"], "_", i)
        if crew in keys(SEQUENCE)
            crew_seq = SEQUENCE[crew]

            # Depart depot
            push!(itemtype, "Start Depot")
            push!(itemidentifier, crew_seq[1])
            push!(optshift, DATE+Minute(SHIFT_START[crew]))
            push!(crewtypeid, crew_row["CrewTypeID"])
            push!(crewtypecode, crew_row["CrewTypeCode"])
            push!(crewnum, i)
            push!(sequence, 1)
            push!(distancetravelforjob, 0)
            push!(pointscalculatedforjob, missing)
            push!(planstartdate, DATE + Minute(SHIFT_START[crew]))
            push!(planenddate, DATE + Minute(SHIFT_START[crew]))

            # Return depot
            push!(itemtype, "Return Depot")
            push!(itemidentifier, crew_seq[end])
            push!(optshift, DATE+Minute(SHIFT_START[crew]))
            push!(crewtypeid, crew_row["CrewTypeID"])
            push!(crewtypecode, crew_row["CrewTypeCode"])
            push!(crewnum, i)
            push!(sequence, length(crew_seq))
            push!(distancetravelforjob, dist_mat[crew_seq[end-1], crew_seq[end]])
            push!(pointscalculatedforjob, missing)
            push!(planstartdate, DATE + Minute(trunc(value(t[crew_seq[end-1]]))) + Second(trunc(60*mod(value(t[crew_seq[end-1]]), 1))) + Minute(trunc(time_mat[crew_seq[end-1], crew_seq[end]])) + Second(trunc(60*mod(time_mat[crew_seq[end-1], crew_seq[end]], 1))))
            push!(planenddate, DATE + Minute(trunc(value(t[crew_seq[end-1]]))) + Second(trunc(60*mod(value(t[crew_seq[end-1]]), 1))) + Minute(trunc(time_mat[crew_seq[end-1], crew_seq[end]])) + Second(trunc(60*mod(time_mat[crew_seq[end-1], crew_seq[end]], 1))))
        end
    end
end

output = DataFrame("ItemType" => itemtype, "ItemIdentifier" => itemidentifier, "OptShift" => optshift, "CrewTypeID" => crewtypeid, "CrewTypeCode" => crewtypecode, "CrewNum" => crewnum, "Sequence" => sequence, "DistanceTravelForJob" => distancetravelforjob, "PointsCalculatedForJob" => pointscalculatedforjob, "PlanStartDate" => planstartdate, "PlanEndDate" => planenddate)

println("OptRunOutput table created successfully.")

##### SAVING SOLUTION #####
# Creating folder structure
if isdir(string("../runs/", request_ID,"/results/X_matrices"))
    rm(string("../runs/", request_ID,"/results/X_matrices"), recursive=true)
    mkpath(string("../runs/", request_ID,"/results/X_matrices"))
else
    mkpath(string("../runs/", request_ID,"/results/X_matrices"))
end

if isdir(string("../runs/", request_ID,"/results/y_matrices"))
    rm(string("../runs/", request_ID,"/results/y_matrices"), recursive=true)
    mkpath(string("../runs/", request_ID,"/results/y_matrices"))
else
    mkpath(string("../runs/", request_ID,"/results/y_matrices"))
end

if isdir(string("../runs/", request_ID,"/results/variables"))
    rm(string("../runs/", request_ID,"/results/variables"), recursive=true)
    mkpath(string("../runs/", request_ID,"/results/variables"))
else
    mkpath(string("../runs/", request_ID,"/results/variables"))
end

# saving
t_df = DataFrame(zip(repeat([request_ID], length(t.axes[1])), t.axes[1], value.(t.data[:])))
rename!(t_df, ["OptRunID", "Identifier", "Time"])
println("t.csv created successfully.")

for k in CREWS
    CSV.write(string("../runs/", request_ID, "/results/X_matrices/X_", k,".csv"), Tables.table(value.(x[:, :, k].data)), header=x.axes[1])
    CSV.write(string("../runs/", request_ID, "/results/y_matrices/y_", k, ".csv"), Tables.table(value.(y.data[:])), header=y.axes[1])
end
CSV.write(string("../runs/", request_ID, "/results/variables/all_variables.csv"), all_var_dict, header=["Variables", "Value"])
CSV.write(string("../runs/", request_ID, "/results/t.csv"), t_df)
CSV.write(string("../runs/", request_ID, "/results/origins.csv"), Tables.table(x.axes[1]))
CSV.write(string("../runs/", request_ID, "/results/destinations.csv"), Tables.table(x.axes[2]))
CSV.write(string("../runs/", request_ID, "/results/sequence.csv"), collect(SEQUENCE), header=["CrewTypeID", "Sequence"])
CSV.write(string("../runs/", request_ID, "/results/optrunmetricsoutput.csv"), DataFrame(metrics), header=true)
CSV.write(string("../runs/", request_ID, "/results/optrunoutput.csv"), output)

println("All results files created and saved successfully.")