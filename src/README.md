# Overview of the source code folder
The source code can be ran either locally on a Downer computer, or on the Roads Data Domain Batch node. There are two execution pathways depending on whether it is a warm start or a normal optimisation run. Warm starts have no defect jobs and only one crew, and is trying to optimise an inspection run. 

 ## execute.py/execute_ws.py
 This is the orchestrator, and is the file that is ran via the terminal. To run this script, it requires two command line options:

 * contract_ID (For example, VRMC_OPS)
 * request_ID (For example, EF4DD65B-B09B-47DB-ABF4-2EB8D9B2A211)
 
 These option parameters are automatically passed down from ADF on Batch. If ran locally, we need to manually specify them. These parameters are then passed down to the subprocesses.

 ~~~
 python execute.py --contract_ID [CONTRACT_ID] --request_ID [REQUEST_ID] 
 ~~~

This scripts does several things:

1. Determines whether it is ran locally or on the Batch node in the load_env_var() function. If it is run locally, all environment variables will point towards the DEV resources on Azure. **The input data is not automatically populated into the DEV storage accounts. So if a run ID is ran on the Batch, and you want to run it locally, you must transfer the data manually into the dmroads/optimisation/input directory.** 
 
2. The extract() function then extracts all the input data based on the environment variables, and saves it locally under the runs/request_ID/data folder. Note that the optimisation start time is currently defaulted to 12am. 
    - v_GetOptimiseInputItems.parquet
    - v_GetOptimiseInspectionInputItems.parquet
    - v_DomainValueDepot.parquet
    - v_GetOptimiseInputCrews.parquet
    - v_GetDistTimeCacheForOptRun.parquet

**Note that we are currently not extracting v_GetOptimiseFunctionDetails.parquet, as it is hard coded in points_curve.py.**

3. After this, it sequentially starts the subprocesses. 

4. If the script is ran on Batch, it will then upload the terminal outputs/errors onto Datalake Storage.

## preprocessing.py/preprocessing_ws.py
This script performs data preprocessing and takes in 3 command line arguments. It does three main things: 
1. Processing the input files into a cleaner format (Possibly redundant)
2. Calculate the score for each job based on the due date. (Based on the start date time command line arguments passed in)
3. Plotting the optimisation run

All files in the preprocessing step gets saved locally in the runs/request_ID/outputs folder.

**Currently the points curve is hard-coded. In the future, v_GetOptimiseFunctionDetails.parquet can be ingested to calculate the score.** However, this has been a deliberate choice. Changing the scoring function can drastically affect model performance. A maximum score of 250 has produced fairly quick and high quality solutions empirically. This is because the score is relative to the travel time (In minutes). A higher maximum score would emphasize completing jobs, while ignoring efficient travel. I have also observed from many runs that a majority of jobs have a minimum score of 50, and not many are in the upper ends of 250.

## optimise.jl/optimise_ws.jl
This is the Julia script that creates the mixed integer linear program model. It ingests the data, and transforms the tables into matrices. The mathematical formulation is under /docs. We are using the JuMP library, and we solve the model using Gurobi.

Once the solve completes, we do some post-processing (On the matrices X, y and t) to produce the results. 
- If it is a warm start, the main result is the t.csv file for the inspection sequence
- If it is a normal optimisation run, the main result is the optrunmetricsoutput.csv and optrunoutput.csv

### Solving method for optimise.jl
For a normal optimisation request, we will normally ingest warms starts for all crews with inspections. This is done via v_GetOptimiseInspectionInputItems. We process this input file as a warm start, but due to rounding errors and how it is stored in the database, ingesting this warm start will produce an infeasible solution (It will violate the constraints by small margin of 0.0001). But the solver can still use this input and re-adjust it slightly to produce a feasible solution. That is why we first do an initial solve with solver parameters SolutionLimit = 1 and MIPFocus = 1. It will normally take the solver less two minutes to re-adjust the values to make a feasible solution.

After it has re-adjusted and produced a feasible solution, we kick off another solve with that solution using the Gurobi NoRel heuristic. NoRel is extremely powerful on improving an existing solution, but is much weaker on finding an initial solution, hence why we need an initial solution first. It will spend a maximum of 1020 seconds in NoRel before stopping.

- Phase 1: MIPFocus = 1 solve to ingest the warm start.
- Phase 2: NoRel heuristic to build from the warm start solution. (1020 seconds)

### Solving method for optimise_ws.jl
We have found that all inspection routes without ramps can be solved. We first do MIPFocus = 1 solving. This will either help find a highly quality solution, or solve to optimiality. If it isn't solved to optimality, we then kick off NoRel heuristic to improve on that solution.

- Phase 1: MipFocus = 1 solve to find a solution and hopefully solve to optimality. (330 seconds)
- Phase 2: NoRel herustic if Phase 1 solve is not optimal. (810 seconds)

## plot_solution.py/plot_ws.py
This plots the results using the folium library. The output is either folium_solution.html for normal optimisations or folium_ws.html for warm starts.

## load.py/load_ws.py
This uploads all the relevant files and results to Datalake Storage. Again, where it gets uploaded to depends on the environment variables configured on the device.
