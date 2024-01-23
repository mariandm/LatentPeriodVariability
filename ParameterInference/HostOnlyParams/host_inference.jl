# Code by Marian Dominguez-Mirazo and David Demory, 2023
# Usage: julia host_inference.jl dataset_folder dataset_id 

# This script infers host-only parameters (mu, growth rate and K, carrying capacity), using host-only growth dynamics.
# It runs 2 steps: Step 1) segments teh data to get a point estimate of mu and K based on growth slope and y axis intersection.
# Step 2) Performs a MCMC using priors informed by step 1).
# The script returns two files: 1. ll_id.csv, and hostchain_id.csv, where id identifies the dynamics data. 

# Read dataset folder and id
tmp_dataid = string(ARGS[1])
tmp_id = string(ARGS[2])

# Load required packages and functions
using DifferentialEquations, CSV, DataFrames, GLM, SparseGrids, Distributions, Dierckx, SpecialFunctions, DelimitedFiles, TimerOutputs, DSP
include("../Functions/segment_estimate.jl")

# Create a folder to store output
tmp_path = tmp_dataid;
mkpath(tmp_path);
# Read dynamics data
dir = "../CreateData/"*tmp_path;
filedir = dir*"/novir_noise_"*tmp_id*".csv"
dataNoVir = DataFrame(CSV.File(filedir,header=false));
# Save timepoints
tsolH = dataNoVir[:,1];
# Get the mean dynamics of three replicates
ysolH = mean(Array(dataNoVir[:,2:4]),dims=2)[:,1];

# Run segment algorithm and save point estimates as csv
println("Running: Segment algorithm");
pars_estimate_segment=segment_estimate(tsolH,ysolH);
writedlm(tmp_path*"/prior_estimates_"*string(tmp_id)*".csv",pars_estimate_segment,",");


###### Running MCMC portion

println("Loading dependencies...")

using MCMCChains, Distributed; 
# Number of cores to be used
addprocs(4);
@everywhere using Turing, Distributions, Statistics, CSV, DataFrames, JLD, DelimitedFiles
@everywhere using LinearAlgebra, DifferentialEquations, ParallelDataTransfer
@everywhere include("../Functions/ODE_SEnIV.jl")

# Allocate space for parallel running
for i in procs()
    @fetch @spawnat i Random.seed!(abs(rand(Int)))
end

# Make relevant information available at every core
sendto(workers(), tmp_id=tmp_id)
sendto(workers(), tmp_path=tmp_path)
sendto(workers(), pars_estimate_segment=pars_estimate_segment)
sendto(workers(), dataNoVir=dataNoVir)
@everywhere id = tmp_id
@everywhere path = tmp_path
@everywhere pars_estimate = pars_estimate_segment

# Get initial conditions
@everywhere filedir = "../CreateData/"*path*"/x0OG_"*id*".csv"
@everywhere x0 = DataFrame(CSV.File(filedir,header=false))


println("Loading dependencies: Done")


# Make previous estimates available in every core
@everywhere begin
    μ_estimate = pars_estimate[1,1]
    K_estimate = pars_estimate[2,1]

    # Create lower bound for the growth rate point estimate
    if μ_estimate < 0.0001
	μ_estimate = 0.01;
    end
end

# Functions to write lognormal in terms of its mode and mean
@everywhere for_mu(mode,mean) = log(mode) + (2/3)*log(mean/mode)
@everywhere for_sigma(mode,mean) = sqrt((2/3)*log(mean/mode))


# Create prior distributions based on point estimates
@everywhere d_μ = truncated(LogNormal(for_mu(μ_estimate,μ_estimate*2),for_sigma(μ_estimate,μ_estimate*2)),0,2)
@everywhere d_K = truncated(LogNormal(for_mu(K_estimate,K_estimate*2),for_sigma(K_estimate,K_estimate*2)),1e3,1e10)

# Set up the problem 
@everywhere begin 
    # Parameters, viral parameters are not relevant for this step of inference
    p = [μ_estimate, K_estimate, 1, 1, 1, 0]
    
    # Initial conditions
    u0=zeros(3,1);
    u0[1] = x0[1,1]
    t0 = 0
    tf = dataNoVir[end,1]
    dt = 0.1
    tmod = collect(t0:dt:tf)
    tspan = (t0,tf)

    # Solve system
    prob = ODEProblem(ODE_SEnIV,u0,tspan,p)
    ymod = solve(prob, Tsit5(),saveat=dt);
    ymod_priors = ymod
end

# Define the MCMC function
@everywhere @model function fitlv(data, prob, dt)
    # Define prior distributions
    μ ~ d_μ
    K ~ d_K
    ϕ = 1
    β = 1
    η = 1
    n = 0
    σ ~ InverseGamma(2,3) #Error likelihood function
    
    # Parameter vector
    pMCMC = [μ, K, ϕ, β, η, n]
    # Solve system
    probMCMC = remake(prob,p=pMCMC,u0=u0)
    ymodMCMC = solve(probMCMC,Tsit5(),saveat=dt)
    # Get indexes where the data timepoints and the simulation match
    idx = indexin(dataNoVir[:,1],ymodMCMC.t)
    # Concatenate dynamics 3 times to compare against 3 replicates
    predicted = vcat(ymodMCMC[1,:][idx],ymodMCMC[1,:][idx],ymodMCMC[1,:][idx])
    
    # Likelihood
    for i =1:length(predicted)
        tmp = predicted[i]
        if tmp<0
            tmp = 0
        end
        log_predicted = log.(tmp)
        data[i,:] ~ Normal(log_predicted, σ)
    end
end

# Define fitting function
@everywhere modelfit = fitlv(log.(vcat(dataNoVir[:,2],dataNoVir[:,3],dataNoVir[:,4])),prob,0.01);

println("Running Host MCMC...")

# Warm-up steps
nburn = 500;
# Steps to save
nstep = 1000;
# Independent chains
nchain = 4;
#Run MCMC
chainfit = reduce(chainscat,pmap(x->sample(modelfit,NUTS(nburn,0.65),nstep,save_state=false),1:nchain))

# Save chain in csv file 
chainarray=Array(chainfit);
writedlm(path*"/hostchain_"*id*".csv",chainarray,",");

