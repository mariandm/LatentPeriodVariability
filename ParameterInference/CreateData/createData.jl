# Code by Marian Dominguez-Mirazo, 2023
# Usage: julia createData.jl dataid

# This script simulates host growth in the absence of virus and 'multi-cycle response curves' for multiple Coefficient of Variation (CV) values, for different parameter datasets. 
# It generates sampled replicates with added noise and saves all data in csv files. 

# Include packages and elevant functions
using DifferentialEquations, Dierckx, Plots, Distributions, DelimitedFiles, SpecialFunctions,DSP
include("../Functions/ODE_SEnIV.jl")
include("../Functions/run_vir.jl")

# Read dataset to create
dataid = string(ARGS[1])
# Create the folder
mkpath(dataid);

# For parameter sources see main text Table 2

# E.coli and lambda-phage parameter set
if dataid=="data1" 
	μ = 1.2;		# Growth rate, hr^-1
	K = 1e8;		# Carrying capacity, CFU/ml
	ϕ = 1e-8;		# Adsorption rate, ml/(CFUxhr)
	β = 200;		# Burst size
	η = 1; 			# Lysis rate, hr^-1
	S_init = 1e5;  	# Initial host density, CFU/ml
	V_init = 1e3; 	# Initial viral density, PFU/ml
	tf_novir = 10; 	# Simulated time for host-only dynamics 
	tf_vir = 12; 	# Simulated time for viral dynamics
	tf_virh = 15; 	# Simulated time for host dynamics

# P. marinus and P-HM2 parameter set
elseif dataid=="data2"
	μ = 0.035; 		# Growth rate, hr^-1
	K = 3e9; 		# Carrying capacity, CFU/ml
	ϕ = 9.3e-10;  	# Adsorption rate, ml/(CFUxhr)
	β = 40; 		# Burst size
	η = 1/5; 		# Lysis rate, hr^-1
	S_init = 1e8; 	# Initial host density, CFU/ml
	V_init = 1e6; 	# Initial viral density, PFU/ml
	tf_novir = 250; # Simulated time for host-only dynamics 
	tf_vir = 30; 	# Simulated time for viral dynamics
	tf_virh = 50; 	# Simulated time for host dynamics

# E. hux and EhV parameter set
elseif dataid=="data3"
	μ = 0.015;  	# Growth rate, hr^-1
	K = 1e9; 		# Carrying capacity, CFU/ml
	ϕ = 1.5e-7;  	# Adsorption rate, ml/(CFUxhr)
	β = 800;  		# Burst size
	η = 1/6;  		# Lysis rate, hr^-1
	S_init = 1e5; 	# Initial host density, CFU/ml
	V_init = 1e3; 	# Initial viral density, PFU/ml
	tf_novir = 1000 # Simulated time for host-only dynamics 
	tf_vir = 30; 	# Simulated time for viral dynamics
	tf_virh = 45; 	# Simulated time for host dynamics
end

CVs = 0.5:-0.05:0.05 		# CV values used to create data 
dt = 0.01 					# Time step
n_points = 10; 				# Number of sampling points
nrep = 3; 					# Number of replicates
noise = 0.3;				# Noise percentage 

# Simulate host growth control without virus present
params = [μ,K,ϕ,β,η,0]; 	
# Initial conditions
u0=zeros(3,1);
u0[1] = S_init; # Susceptible host
# Simulate dynamics
tspan = (0,tf_novir)
prob = ODEProblem(ODE_SEnIV,u0,tspan,params);
y_novir = solve(prob, Tsit5(),saveat=dt);
ymod_novir = hcat(y_novir.t,reverse(rotl90(hcat(y_novir.u...)),dims=1));


# Run simulations for each CV value
for i = 1:size(CVs)[1]

	# Update parameters
	cv = CVs[i];
	p = [μ,K,ϕ,β,η,cv]; 

	# Run simulation
	@time ymod_vir = run_vir(p,[S_init,V_init],maximum([tf_vir,tf_virh]),dt);
	
	# Sample and add noise
	
	# Sample host growth control
	# Use timepoints that are half hours or full hours (to simulate real experimental sampling)
	half_hours = findall(x->mod(x,0.5)==0,ymod_novir[:,1]);
	# Get the indexes of n_points time points equally distant 
	idx = half_hours[Int64.(round.(vcat(1, (1:n_points) .* (size(half_hours)[1] - 1)/n_points .+ 1)))];
	# Save time points and host density of the respective indexes
	novir_sample = hcat(ymod_novir[idx,1],ymod_novir[idx,2]);

	# Sample viral dynamics
	# Get timepoints that are half hours
	half_hours = findall(x->mod(x,0.5)==0,ymod_vir[ymod_vir[:,1] .<= tf_vir,1]);
	# Get indexes of equally distant half hours
	idx = half_hours[Int64.(round.(vcat(1, (1:n_points) .* (size(half_hours)[1] - 1)/n_points .+ 1)))];
	# Save time points and viral densities
	vir_samplev = hcat(ymod_vir[idx,1],ymod_vir[idx,end]);

	# Sample total host dynamics
	# Get timepoints where total host density is higher than 1
	above1 = sum(ymod_vir[:,2:end-1],dims=2).>1 .* ymod_vir[:,1] .<= tf_virh;
	# Get timepoints that are haf hours
	half_hours = findall(x->mod(x,0.5)==0,ymod_vir[above1,1]);
	# Get indexes of equally distant half hours
	idx = half_hours[Int64.(round.(vcat(1, (1:n_points) .* (size(half_hours)[1] - 1)/n_points .+ 1)))];
	# Save time points and total host densities
	vir_sampleh = hcat(ymod_vir[idx,1],sum(ymod_vir[idx,2:end-1],dims=2));

	# Create replicates with added noise

	# Save timepoints in the first column
	novir_noise_array = novir_sample[:,1];
	vir_noisev_array = vir_samplev[:,1];
	vir_noiseh_array = vir_sampleh[:,1];

	# Create nrep replicates with added noise
	# Save replicates in the following columns
	for r = 1:nrep
		novir_noise_array = hcat(novir_noise_array, noise .* rand(Normal(),size(novir_sample)[1]) .* novir_sample[:,2] .+ novir_sample[:,2]);
		vir_noisev_array = hcat(vir_noisev_array, noise .* rand(Normal(),size(vir_samplev)[1]) .* vir_samplev[:,2] .+ vir_samplev[:,2]);
		vir_noiseh_array = hcat(vir_noiseh_array, noise .* rand(Normal(),size(vir_sampleh)[1]) .* vir_sampleh[:,2] .+ vir_sampleh[:,2]);
	end

	# Save sampled replicates as csv files
	# Host growth control
	writedlm(dataid*"/novir_noise_"*string(i)*".csv",novir_noise_array,",");
	# Viral dynamics
	writedlm(dataid*"/vir_noisev_"*string(i)*".csv",vir_noisev_array,",");
	# Host dynamics
	writedlm(dataid*"/vir_noiseh_"*string(i)*".csv",vir_noiseh_array,",");


	# Save original values and parameters
	# Initial conditions
	writedlm(dataid*"/x0OG_"*string(i)*".csv",[S_init,V_init],",");
	# Parameters
	writedlm(dataid*"/parsOG_"*string(i)*".csv",[μ,K,ϕ,β,η,cv],",");
	# Host growth control full dynamics 
	writedlm(dataid*"/novirOG_"*string(i)*".csv",ymod_novir,",");
	# Host and viral dynamics
	writedlm(dataid*"/virOG_"*string(i)*".csv",ymod_vir,",");

end
