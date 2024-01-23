# Code by Marian Dominguez-Mirazo, 2023
# Usage: julia step1_likelihood.jl dataset_folder dataset_id 

# This script is the first step for virus-host parameter inference. 
# It performs an iterative grid search in a range of biologically-relevant parameters and returns tow files.
# The first file contains the maximum likelihood point estimates and second one, the parameter space with the likelihood associated values. 

# Read datset folder and id
dataid = string(ARGS[1])
id = string(ARGS[2])

# Load dependencies
using DifferentialEquations, CSV, DataFrames, GLM, SparseGrids, Distributions, Dierckx, SpecialFunctions, DelimitedFiles, TimerOutputs, DSP
include("../Functions/prior_estimate_ll.jl")


dt = 0.01;
Tol = 0.01;
cvTol = 0.049;

# Read previously estimated host parameters (see ../HostOnlyParams/host_inference.jl)
filedir = "../HostOnlyParams/"*dataid*"/hostchain_"*id*".csv";
host_chains = DataFrame(CSV.File(filedir,header=false));
pars_estimate_host = [mean(host_chains[:,1]),mean(host_chains[:,2])];

# Read Initial conditions
dir = "../CreateData/"*dataid
filedir = dir*"/x0OG_"*id*".csv"
x0 = DataFrame(CSV.File(filedir,header=false))[:,1];

# Read replicates data
# Viral dynamics
filedir = dir*"/vir_noisev_"*id*".csv"
rep_VirV = DataFrame(CSV.File(filedir,header=false));
# Host dynamics
filedir = dir*"/vir_noiseh_"*id*".csv"
rep_VirH = DataFrame(CSV.File(filedir,header=false));

# Calculate standard deviation of the three replicates
std_VirV = std.(collect(eachrow(reduce(hcat, [log.(rep_VirV[:,2]),log.(rep_VirV[:,3]),log.(rep_VirV[:,4])]))));
std_VirH = std.(collect(eachrow(reduce(hcat, [log.(rep_VirH[:,2]),log.(rep_VirH[:,3]),log.(rep_VirH[:,4])]))));

# Setup structure to consider three replicates
# Viral dynamics standard deviation
std_VirV = vcat(std_VirV,std_VirV,std_VirV); 
# Total host standard deviation
std_VirH = vcat(std_VirH,std_VirH,std_VirH); 
# Time points for total host
tsolHV = vcat(rep_VirH[:,1],rep_VirH[:,1],rep_VirH[:,1]); 
# Total host dynamics
ysolHV = vcat(rep_VirH[:,2],rep_VirH[:,3],rep_VirH[:,4]); 
# Time points for viral counts
tsolV = vcat(rep_VirV[:,1],rep_VirV[:,1],rep_VirV[:,1]);
# Viral dynamics
ysolV = vcat(rep_VirV[:,2],rep_VirV[:,3],rep_VirV[:,4]);
    
    
# The code can estimate all parameters, or assume that some parameters are known beforehand. 
# The main text estimates all parameters
# To assume some parameters are known, uncomment the line that best describes the desired situation
prevparam = [NaN,NaN]; 						# Estimate all parameters

# Get parameters with which data was created
#filedir = dir*"/parsOG_"*id*".csv"
#parsOG = DataFrame(CSV.File(filedir,header=false))[:,1];
# Inform the program which parameters are known
#prevparam = [NaN,parsOG[4]*1];	  			# Assume burst size is known, estimate phi, eta, cv
#prevparam = [parsOG[3]*1,NaN];	            # Assume phi is known, estimate beta, eta, cv
#prevparam = [parsOG[3]*1,parsOG[4]*1];    	# Assume phi and beta are known, estimate eta and cv

# Run parameter estimation
pars_estimate, ll_map = prior_estimate(pars_estimate_host,tsolHV,ysolHV,tsolV,ysolV,std_VirV,std_VirH,x0,prevparam,dt);

# Create output folder
path="step1_likelihood/"*dataid;
mkpath(path);
# Save point estimates
writedlm(path*"/estimates_"*string(id)*".csv",pars_estimate,",");
# Save likelihood neighbor
writedlm(path*"/ll_"*string(id)*".csv",  ll_map, ',');
