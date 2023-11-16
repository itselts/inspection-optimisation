# Automation and optimisation of road network maintenance
An automated decision support tool that optimizes the scheduling and routing of road network maintenance. It takes in both routine inspections as well as point defects. The underlying model is a mixed integer linear program created in Julia's JuMP library, and is solved using Gurobi. The model formulation is in the /docs directory, and can be classified as a capacitated, windy, multi-depot, genenral routing problem. Developed by Elton Shi and David Ming, the first version was released in February 2023.

![alt text](https://github.com/itselts/inspection-optimisation/tree/main/runs/7B81F1F2-2BB6-45EA-BA70-DD823385BD0A/results)

## Model features
- Multi-depot
- Non-homogenous crews
    - Crew job constraints
    - Different start and end depots
    - Routine inspections crews and/or point defect crews
- Due dates
    - Penalties
- Road directionality

## How to run
There are numerous things that must be setup first before one is able to run an optimisation. This repository is integrated and setup specifically to be run on the Downer network (AzureDefaultCredential), either locally, or directly on Azure Batch. A .env file is required with the relevant environment variables pointing to the correct Azure resource URLs (Specifically Batch Storage Account). 

### Clone the repository
Ensure you have Git installed on your computer. For this guide, we will specifically use Github and VSCode. (https://learn.microsoft.com/en-us/azure/developer/javascript/how-to/with-visual-studio-code/clone-github-repository?tabs=create-repo-command-palette%2Cinitialize-repo-activity-bar%2Ccreate-branch-command-palette%2Ccommit-changes-command-palette%2Cpush-command-palette)

1. Install Github Pull Requests and Issues extension in VSCode.
2. Open the command palette and Git Clone https://github.com/DownerRoads/routine_inspection_optimisation_project.git
3. Save the remote repository in a suitable location.

### Dependencies
#### Python
Ensure you have Python installed on your computer. It is recommended to create a Python virtual environment in the root folder. There is no preference to how you create your own virtual environment. In this guide we will use the venv module in Python at https://docs.python.org/3/library/venv.html. 

1. Navigate to the root folder of the repository in a terminal window.
2. To create a virtual environment, enter in the terminal: 
~~~
python -m venv .venv
~~~
3. To activate the virtual environment, enter in the terminal: 
~~~
.venv/Scripts/activate.bat 
~~~
4. To install the packages listed in requirements.txt, enter in the terminal: 
~~~
pip install -r requirements.txt
~~~

#### Julia
Download Julia at https://julialang.org/downloads/. Julia dependencies are in the Manifest.toml and Project.toml files in the /src directory. The following guides explains how to install Julia dependencies locally.
- https://pkgdocs.julialang.org/v1/toml-files/ 
- https://stackoverflow.com/questions/72061475/installing-julia-packages-using-a-toml-file 

#### Gurobi
This tool uses Gurobi as its solver, hence an installation of Gurobi and a valid licence is required. Login and download from https://www.gurobi.com/downloads/gurobi-software/. The Gurobi library in Julia (https://juliapackages.com/p/gurobi) then acts as a wrapper when used in conjunction with the JuMP library, hence nothing else needs to be done. If we were using Python, we would need to use the special gurobipy library and its specific modelling syntax.

### Executing an optimisation run
To run an optimisation for a specified request ID and contract ID:

1. Activate the correct Python virtual environment.
2. Set the active directory to src
3. Run execute.py/execute_ws.py script via the command line with the following command:
~~~
python execute.py --contract_ID [CONTRACT_ID] --request_ID [REQUEST_ID]
~~~

If ran locally, it will attempt to access:
- Storage account: https://adlsdevrdsddaze02.dfs.core.windows.net
- Container: data
- Directory: dmroads/optimisation/input/CONTRACT_ID/REQUEST_ID

#### Normal optimisation example
~~~
python execute.py --contract_ID VRMC_OPS --request_ID 7B81F1F2-2BB6-45EA-BA70-DD823385BD0A
~~~
- 203 inspection segments
- 762 jobs
- 2 depots
- 2 crews

#### Warm start optimisation example
~~~
python execute_ws.py --contract_ID VRMC_OPS --request_ID B6EDB008-8717-4988-8AD9-DE9C8BB43B2A
~~~
- Wodonga depot with 105 inspection segments

## Future improvements (From easiest to hardest)
- Improve the scoring function. Currently, a majority of jobs have a score of 50. Possibly need a steeper gradient too.
- Ingest v_GetOptimiseFunctionDetails. Remove hard coding.
- Job overload. There needs to be some vetting of jobs, either through date filtering or geofencing. If a crew can only do about 15 jobs a day, we don't need to optimise for 500+ jobs. Resource intensive and drastically reduces solution quality/speed.
- Intake crew travel cost into the model.
- Introduce penalties for having a job go overdue.
- Continuing a solve to improve on a current solution.
- Absorb jobs on inspection roads into the inspection edge. 
    - Major alteration to the model, but fixes the biggest deficiency of the current model. By doing this, we remove "backtracking" from exiting an inspection, then re-entering the inspection. The model currently might not even have the ability to do jobs on freeways as there is no way to "backtrack". 
    - Can also remove the need to segment inspections into such high granularity/fidelity. Drastically reduce model size.
- Multi-day horizon optimisations. Could just run multiple optimisations in series.
- Research into creating an internal server to calculate the travel matrices data instead of querying Azure.
