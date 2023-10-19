This folder was the first proof of concept developed by David Ming to determine the feasibility of optimising for point defects and inspection segments together. 

## make_network.jl 
A script that creates an interactive window for the user to make their own dummy network problem. It allows the user to create inspection roads as well as putting in point defect jobs.

## linear_task.jl
A script that solves the network created from make_network.jl

## test_network_x.json
User created networks from make_network.jl. It is formatted such that it can be ingested by linear_task.jl.

## Picture1.png
This shows the output of solving one of the dummy example optimally. 