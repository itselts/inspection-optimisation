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

# request_ID = "22E82CDB-420D-4D98-ACE5-AA59728F59FE"
# s_time_str = "00:00"
# s_date_str = "2023-06-16"
# DATE = DateTime(s_date_str)
# DOW = dayname(DATE)
# TIME = Time(s_time_str)

###### DATA INGESTION #####
print("Ingesting and transforming data...")

# inspection segments
seg_df = @pipe read_parquet(string("../runs/", request_ID , "/data/v_GetOptimiseInspectionInputItems.parquet"))  |> DataFrame

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

pairwise1 = copy(pairwise_cache[!, [:PointAID_FK, :PointBID_FK, :Distance, :Time]])
pairwise1.PointAID_FK = string.(pairwise1.PointAID_FK) # Type conversions
pairwise1.PointBID_FK = string.(pairwise1.PointBID_FK)
pairwise1.Time = Float64.(pairwise1.Time)
pairwise1.Distance = Float64.(pairwise1.Distance)

pairwise2 = seg_df[!, [:SegStartIdentifier, :SegEndIdentifier, :EstQty, :EstDuration]] # segments pairwise
pairwise2.EstDuration = Float64.(pairwise2.EstDuration) # Type conversions
pairwise2.EstQty = Float64.(pairwise2.EstQty) ./ 1000 # Type conversion, then to kilometres
rename!(pairwise2, :SegStartIdentifier => :PointAID_FK, :SegEndIdentifier => :PointBID_FK, :EstQty => :Distance, :EstDuration => :Time)

pairwise = vcat(pairwise1, pairwise2) # Full pairwise list 

dist_df = @pipe unstack(pairwise, :PointAID_FK, :PointBID_FK, :Distance) # Pairwise to DataFrame 
time_df = @pipe unstack(pairwise, :PointAID_FK, :PointBID_FK, :Time)

# NODES FROM THE CACHE TABLE (c.f NODES FROM THE INPUT ITEM TABLES)
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
    if crew_num == 0
        println("********** There are no crews rostered on ", DOW, " for ", row["CrewTypeCode"] ,". **********")
    end
    for i in 1:crew_num
        crew = string(row["CrewTypeCode"], "_", i)
        start_time = Time(row[string(DOW, "StartTime")] - Dates.Hour(12))
        end_time = Time(row[string(DOW, "EndTime")] - Dates.Hour(12))
        SHIFT_START[crew] = (60*Hour(start_time).value + Minute(start_time).value)
        SHIFT_END[crew] = (60*Hour(end_time).value + Minute(end_time).value) 
    end
end

##### SETS FROM INPUT DATA #####
NODES_DATA = Dict()
for row in eachrow(seg_df)
    NODES_DATA[string(row["SegStartIdentifier"])] = (parse(Float64, row["LongStart"]), parse(Float64, row["LatStart"]))
    NODES_DATA[string(row["SegEndIdentifier"])] = (parse(Float64, row["LongEnd"]), parse(Float64, row["LatEnd"]))
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

JOBS = SEGMENT_POINTS
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
    end
end


##### CONSTANTS/MODEL INPUTS #####
REVERSIBLE = Dict(string(row["ItemIdentifier"]) => row["Reversible"] for row in eachrow(seg_df)) # Dictionary comprehension

println("ok")

##### MODELLING #####
model = Model(Gurobi.Optimizer)

print("Creating the model...")

@variable(model, x[NODES, NODES, CREWS], Bin)
@variable(model, y[JOBS, CREWS], Bin)
@variable(model, t[JOBS] >= 0)

# Big M constants
M1 = maximum(values(SHIFT_END)) + maximum(skipmissing(time_mat.data))
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
    
    
    # @constraint(model, sum(x[i, j, k] for k in CREWS) == 1) # Enforcing forward travel direction only for warm start
    # @constraint(model, sum(x[j, i, k] for k in CREWS) == 0)
    
end

# 2a. conservation of flow. If a crew enters, crew must also leave 
for j in NODES
    for k in CREWS
        @constraint(model, sum(x[i, j, k] for i in NODES) == sum(x[j, i, k] for i in NODES))
    end
end

# 2b. at most one crew can be active leaving jobs
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
    try  # Check if entry in time matrix is missing
        @constraint(model, t[j] >= t[i] + time_mat[i,j] - M1*(1 - sum(x[i,j,k] for k in CREWS))) 
    catch ex
        continue
    end
end

# 10. shift start and end time
for j in JOBS
    try # Check if entry in time matrix is missing
        begin @constraint(model, t[j] <= sum(x[j, crew_start_end[k]["end"], k]*(SHIFT_END[k] - time_mat[j, crew_start_end[k]["end"]]) for k in CREWS) 
            + M2*(1 - sum(x[j, crew_start_end[k]["end"], k] for k in CREWS))) end # End time
    catch ex
    end
    try
        @constraint(model, t[j] >= sum(x[crew_start_end[k]["start"], j, k]*(SHIFT_START[k] + time_mat[crew_start_end[k]["start"], j]) for k in CREWS)) # Start time
    catch ex
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
@expression(model, tot_time, sum(coalesce(time_mat[i, j], 5000) * x[i, j, k] for i in NODES, j in NODES, k in CREWS))
@expression(model, tot_jobs, sum(y[i, k] for i in JOBS, k in CREWS))
@objective(model, Min, 100*tot_time)

println("ok")

##### SOLVING #####
set_optimizer_attribute(model, "TimeLimit", 330)
#set_optimizer_attribute(model, "SolutionLimit", 1)
set_optimizer_attribute(model, "MIPFocus", 1)

optimize!(model)

##### SAVING INITIAL SOLVE ######
# variable values
all_var_dict = Dict(zip(all_variables(model), value.(all_variables(model))))

if isdir(string("../runs/", request_ID,"/warm_start/variables"))
    rm(string("../runs/", request_ID,"/warm_start/variables"), recursive=true)
    mkpath(string("../runs/", request_ID,"/warm_start/variables"))
else
    mkpath(string("../runs/", request_ID,"/warm_start/variables"))
end
CSV.write(string("../runs/", request_ID, "/warm_start/variables/all_variables.csv"), all_var_dict, header=["Variable", "Value"]) # Saving solve results.


##### WARM START WITH A PRIOR SOLVE #####
if relative_gap(model) >= 0.00001
    println("")
    println("Starting phase 2 solve to improve on the solution...")
    println("")
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
            set_start_value(t[i], val) 
        end
    end

    MOIU.reset_optimizer(model)
    
    set_optimizer_attribute(model, "SolutionLimit", 200) # Can't reset optimizer attribute to null
    set_optimizer_attribute(model, "TimeLimit", 810) # Can't reset optimizer attribute to null
    set_optimizer_attribute(model, "NoRelHeurTime", 810)
    
    optimize!(model)
end


# variable values
all_var_dict = Dict(zip(all_variables(model), value.(all_variables(model))))
CSV.write(string("../runs/", request_ID, "/warm_start/variables/all_variables.csv"), all_var_dict, header=["Variable", "Value"]) # Updating solve results.

##### PROCESSING ORDER OF POINTS ######
SEQUENCE = Dict{String, Vector{String}}()

for k in CREWS
    i = 1
    start = crew_start_end[k]["start"]
    if !(isnothing(findfirst([round(value(x[start, j, k])) for j in NODES] .!= 0))) # Finding if the crew leaves the start depot
        SEQUENCE[k] = []
        to = @pipe findfirst([round(value(x[start, j, k])) for j in NODES] .!= 0) |> NODES[_] # Find the numeric index of the node that it leaves to, then get the node label.
        push!(SEQUENCE[k], start)
        push!(SEQUENCE[k], to)
        while to != crew_start_end[k]["end"]
            from = to
            to = @pipe findfirst([round(value(x[from, j, k])) for j in NODES] .!= 0) |> NODES[_]
            push!(SEQUENCE[k], to)
        end
    end
end


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
                push!(planstartdate, DATE + Minute(trunc(value(t[start]))) + Second(trunc(60*mod(value(t[start]), 1)))) # RESCALING
                push!(planenddate, DATE + Minute(trunc(value(t[start]))) + Second(trunc(60*mod(value(t[start]), 1)))) # RESCALING
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
        crew_seq = SEQUENCE[crew]

        # Depart depot
        push!(itemtype, "Start Depot")
        push!(itemidentifier, crew_seq[1])
        push!(optshift, DATE+Minute(SHIFT_START[crew])) # RESCALING
        push!(crewtypeid, crew_row["CrewTypeID"])
        push!(crewtypecode, crew_row["CrewTypeCode"])
        push!(crewnum, i)
        push!(sequence, 1)
        push!(distancetravelforjob, 0)
        push!(pointscalculatedforjob, missing)
        push!(planstartdate, DATE + Minute(SHIFT_START[crew])) # RESCALING
        push!(planenddate, DATE + Minute(SHIFT_START[crew])) # RESCALING

        # Return depot
        push!(itemtype, "Return Depot")
        push!(itemidentifier, crew_seq[end])
        push!(optshift, DATE+Minute(SHIFT_START[crew])) # RESCALING
        push!(crewtypeid, crew_row["CrewTypeID"])
        push!(crewtypecode, crew_row["CrewTypeCode"])
        push!(crewnum, i)
        push!(sequence, length(crew_seq))
        push!(distancetravelforjob, dist_mat[crew_seq[end-1], crew_seq[end]])
        push!(pointscalculatedforjob, missing)
        push!(planstartdate, DATE + Minute(trunc(value(t[crew_seq[end-1]]))) + Second(trunc(60*mod(value(t[crew_seq[end-1]]), 1))) + Minute(trunc(time_mat[crew_seq[end-1], crew_seq[end]])) + Second(trunc(60*mod(time_mat[crew_seq[end-1], crew_seq[end]], 1)))) # RESCALING
        push!(planenddate, DATE + Minute(trunc(value(t[crew_seq[end-1]]))) + Second(trunc(60*mod(value(t[crew_seq[end-1]]), 1))) + Minute(trunc(time_mat[crew_seq[end-1], crew_seq[end]])) + Second(trunc(60*mod(time_mat[crew_seq[end-1], crew_seq[end]], 1)))) # RESCALING
    end
end

output = DataFrame("ItemType" => itemtype, "ItemIdentifier" => itemidentifier, "OptShift" => optshift, "CrewTypeID" => crewtypeid, "CrewTypeCode" => crewtypecode, "CrewNum" => crewnum, "Sequence" => sequence, "DistanceTravelForJob" => distancetravelforjob, "PointsCalculatedForJob" => pointscalculatedforjob, "PlanStartDate" => planstartdate, "PlanEndDate" => planenddate)

println("OptRunOutput table created successfully.")

##### SAVING SOLUTION #####
# Creating folder structure
if isdir(string("../runs/", request_ID,"/warm_start/X_matrices"))
    rm(string("../runs/", request_ID,"/warm_start/X_matrices"), recursive=true)
    mkpath(string("../runs/", request_ID,"/warm_start/X_matrices"))
else
    mkpath(string("../runs/", request_ID,"/warm_start/X_matrices"))
end

if isdir(string("../runs/", request_ID,"/warm_start/y_matrices"))
    rm(string("../runs/", request_ID,"/warm_start/y_matrices"), recursive=true)
    mkpath(string("../runs/", request_ID,"/warm_start/y_matrices"))
else
    mkpath(string("../runs/", request_ID,"/warm_start/y_matrices"))
end
    
# saving
t_df = DataFrame(zip(repeat([request_ID], length(t.axes[1])), t.axes[1], (value.(t.data[:]) .- SHIFT_START[CREWS[1]]))) # RESCALING BY 100
rename!(t_df, ["OptRunID", "Identifier", "Time"])

for k in CREWS
    CSV.write(string("../runs/", request_ID, "/warm_start/X_matrices/X_", k,".csv"), Tables.table(value.(x[:, :, k].data)), header=x.axes[1])
    CSV.write(string("../runs/", request_ID, "/warm_start/y_matrices/y_", k, ".csv"), Tables.table(value.(y.data[:])), header=y.axes[1])
end
CSV.write(string("../runs/", request_ID, "/warm_start/t.csv"), t_df)
CSV.write(string("../runs/", request_ID, "/warm_start/origins.csv"), Tables.table(x.axes[1]))
CSV.write(string("../runs/", request_ID, "/warm_start/destinations.csv"), Tables.table(x.axes[2]))
CSV.write(string("../runs/", request_ID, "/warm_start/sequence.csv"), collect(SEQUENCE), header=["CrewTypeID", "Sequence"])
CSV.write(string("../runs/", request_ID, "/warm_start/optrunoutput.csv"), output)

println("All results files created and saved successfully.")