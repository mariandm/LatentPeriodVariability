# Code by Marian Dominguez-Mirazo and David Demory, 2023
# Usage: julia step2_round2_MCMC_vir.jl dataset_folder dataset_id chain_id
# This script takes previously calculated first round of MCMC (see step2 round 1) to inform priors for MCMC parameter estimation.
# It returns a csv file with the chains in the order: phi, beta, eta, cv, sigma_h, sigma_v

# Read dataset folder, dataset id, and chain identifier
tmp_dataid = string(ARGS[1])
tmp_id = string(ARGS[2])
tmp_runid = string(ARGS[3])
println("Loading dependencies...")

# Load dependencies and allocate cores
using MCMCChains, Distributed; addprocs(1)
@everywhere using Turing, Distributions, Statistics, CSV, DataFrames, JLD, DelimitedFiles
@everywhere using LinearAlgebra, DifferentialEquations, ParallelDataTransfer, Dierckx, DSP, SpecialFunctions
@everywhere include("../Functions/ODE_SEnIV.jl")
@everywhere include("../Functions/run_vir.jl")

for i in procs()
    @fetch @spawnat i Random.seed!(abs(rand(Int)))
end

# Make information available in all cores
sendto(workers(), tmp_id=tmp_id)
sendto(workers(), tmp_dataid=tmp_dataid)
sendto(workers(), tmp_runid=tmp_runid)
@everywhere id = tmp_id
@everywhere dataid = tmp_dataid
@everywhere runid = tmp_runid
@everywhere path = "step2_MCMC/round2/"*dataid

# Create storage folder
mkpath(path);

# Load data replicates
@everywhere filedir = "../CreateData/"*dataid*"/vir_noisev_"*id*".csv"
@everywhere dataVir = DataFrame(CSV.File(filedir,header=false))
@everywhere filedir = "../CreateData/"*dataid*"/vir_noiseh_"*id*".csv"
@everywhere dataVirH = DataFrame(CSV.File(filedir,header=false))

# Load initial conditions
@everywhere filedir = "../CreateData/"*dataid*"/x0OG_"*id*".csv"
@everywhere x0 = DataFrame(CSV.File(filedir,header=false))

# Load initial parameter estimates (step 1)
@everywhere filedir = "step1_likelihood/"*dataid*"/estimates_"*id*".csv"
@everywhere pars_estimate = DataFrame(CSV.File(filedir,header=false))

# Load first round of chains (step 2 round 1)
@everywhere filedir = "step2_MCMC/round1/"*dataid*"/viruschain_"*id*"_round1.csv";
@everywhere round1_chains = DataFrame(CSV.File(filedir,header=false));

# Time step
@everywhere dt = 0.1;

println("Loading dependencies: Done")

# Make previous estimates available in every core
@everywhere begin
    μ_estimate2 = pars_estimate[1,1];
    K_estimate2 = pars_estimate[2,1];
    ϕ_estimate = mean(round1_chains[:,1]);
    ϕ_std = std(round1_chains[:,1]);
    β_estimate = mean(round1_chains[:,2]);
    β_std = std(round1_chains[:,2]);
    η_estimate = mean(round1_chains[:,3]);
    η_std = std(round1_chains[:,3]);
    cv_estimate = mean(round1_chains[:,4]);
    cv_std = std(round1_chains[:,4]);

    # Get simulation time
    tf = max(dataVir[end,1],dataVirH[end,1])
end

# Functions to write lognormal in terms of its mode and mean
@everywhere for_mu(mode,mean) = log(mode) + (2/3)*log(mean/mode)
@everywhere for_sigma(mode,mean) = sqrt((2/3)*log(mean/mode))

# Create prior distributions based on first round of MCMC
@everywhere d_ϕ = truncated(Normal(ϕ_estimate,2*ϕ_std),-11,-5)
@everywhere d_β = truncated(LogNormal(for_mu(β_estimate,β_estimate+2*β_std),for_sigma(β_estimate,β_estimate+2*β_std)),1,1000)
@everywhere d_η = truncated(LogNormal(for_mu(η_estimate,η_estimate+2*η_std),for_sigma(η_estimate,η_estimate+2*η_std)),1/12,60/10)
@everywhere d_CV = truncated(LogNormal(for_mu(cv_estimate,cv_estimate+2*cv_std),for_sigma(cv_estimate,cv_estimate+2*cv_std)),0.03,1) 

# Define the MCMC function
@everywhere @model function fitlv(data_h, data_v,tf,dt)
    # Assume previously estimated host-only parameters
    μ = μ_estimate2
    K = K_estimate2
    # Assign prior distributions
    ϕ ~ d_ϕ
    β ~ d_β
    η ~ d_η
    cv ~ d_CV
    
    σ_h ~ truncated(InverseGamma(1,0.5),0,2) #sigma likelihood function for host data
    σ_v ~ truncated(InverseGamma(1,0.5),0,2) #sigma likelihood function for viral data

    # Assign parameter vector
    # Using Turing.ForwardDiff.value allows us to sample cv and transform into the best approximation of integer n 
    p = [μ,K,10^Turing.ForwardDiff.value(ϕ),Turing.ForwardDiff.value(β),Turing.ForwardDiff.value(η),Turing.ForwardDiff.value(cv)]
    
    # Run simulation
    ymod = run_vir(p,[x0[1,1],x0[end,1]],tf,dt);

    # Get indexes where the data timepoints and the simulation match
    idxv = indexin(dataVir[:,1],ymod[:,1])
    idxh = indexin(dataVirH[:,1],ymod[:,1])
    # Concatenate dynamics 3 times to compare against 3 replicates
    predicted_v = vcat(ymod[:,end][idxv], ymod[:,end][idxv], ymod[:,end][idxv]);
    predicted_h = vcat(sum(ymod[:,2:end-1],dims=2)[idxh], sum(ymod[:,2:end-1],dims=2)[idxh], sum(ymod[:,2:end-1],dims=2)[idxh]);

    #likelihood host
    for i =1:length(predicted_h)
        tmp_h = predicted_h[i]
        # Prevent values below 0
        if tmp_h<0
            tmp_h = 0
        end
        log_predicted_h = log.(tmp_h)
        data_h[i,:] ~ Normal(log_predicted_h, σ_h)
    end
    
    #likelihood virus
    for i =1:length(predicted_v)
        tmp_v = predicted_v[i]
        # Prevent values below 0
        if tmp_v<0
            tmp_v = 0
        end
        log_predicted_v = log.(tmp_v)
        data_v[i,:] ~ Normal(log_predicted_v, σ_v)
    end
end

# Define fitting function
@everywhere modelfit = fitlv(log.(vcat(dataVirH[:,2],dataVirH[:,3],dataVirH[:,4])),log.(vcat(dataVir[:,2],dataVir[:,3],dataVir[:,4])),tf,dt)

println("Running Viral MCMC...")

# Warm-up steps
nburn = 1000
# Steps to save
nstep = 1000
# Number of chains
nchain = 1
#Run MCMC
chainfit = reduce(chainscat,pmap(x->sample(modelfit,NUTS(nburn,0.45),nstep,save_state=false),1:nchain))

# Save chain in csv file 
chainarray=Array(chainfit);
writedlm(path*"/viruschain"*id*"_round2_"*runid*".csv",chainarray,",");

println("Running Viral MCMC: Done")
