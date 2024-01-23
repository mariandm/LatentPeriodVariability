# Code by Marian Dominguez-Mirazo, 2023

# This function runs an iterative search for the combination of parameters that reduces the error from data fit with maximum likelihood.

# Load dependencies
include("../Functions/likelihood_matrix.jl")

function prior_estimate(pars_estimate_host,tsolHV,ysolHV,tsolV,ysolV,ste_VirV,ste_VirH,x0,prevparam,dt)
	# pars_estimate_host: Host parameters previously estimated
	# tsolHV: Total host timepoints
	# ysolHV: Total host dynamics
	# tsolV: Viral timepoints
	# ysolV: Viral dynamics
	# ste_VirV: Standard error for viral dynamics
	# ste_VirH: Standard error for total host dynamics
	# x0: Initial conditions
	# prevparams: Previously known parameters
	# dt: Time step

	pars_estimate_1 = pars_estimate_host;
	push!(pars_estimate_1,prevparam[1]); push!(pars_estimate_1,prevparam[2]);

	# Determine based on prevparam input if phi is previously known, if beta is knwon, if both are, or neither
	# See vector structure in ../VirusHostParams/step1_likelihood.jl
	situation = sum(.!isnan.(prevparam) .* [2,1]); #turning binary conditions into decimal

	#Iterative likelihood matrix
	println("Running: Likelihood matrix");
	println("This might take a bit");
	println("Patience is a virtue...");

	# Run three iterations of the matrix search, restricting the search space every step (with neis)
	estimates, neis, chi_values = likelihood_matrix(pars_estimate_1,tsolV,ysolV,tsolHV,ysolHV,ste_VirV,ste_VirH,x0,situation,([],[],[],[]),dt,"log");
	estimates, neis, chi_values = likelihood_matrix(pars_estimate_1,tsolV,ysolV,tsolHV,ysolHV,ste_VirV,ste_VirH,x0,situation,(neis[1],neis[2],neis[3],neis[4]),dt,"linear");
	estimates, neis, chi_values = likelihood_matrix(pars_estimate_1,tsolV,ysolV,tsolHV,ysolHV,ste_VirV,ste_VirH,x0,situation,(neis[1],neis[2],neis[3],neis[4]),dt,"linear");


	if situation==0 		# Estimate all parameters
		pars_estimate = [pars_estimate_1[1],pars_estimate_1[2],estimates[3],estimates[4],estimates[1],estimates[2]]
	elseif situation==1 	# Assume burst size is known, estimate phi, eta, cv
		pars_estimate = [pars_estimate_1[1],pars_estimate_1[2],estimates[3],prevparam[2],estimates[1],estimates[2]]
	elseif situation==2		# Assume phi is known, estimate beta, eta, cv
		pars_estimate = [pars_estimate_1[1],pars_estimate_1[2],prevparam[1],estimates[4],estimates[1],estimates[2]]
	else 					# Assume phi and beta are known, estimate eta and cv
		pars_estimate = [pars_estimate_1[1],pars_estimate_1[2],prevparam[1],prevparam[2],estimates[1],estimates[2]]
	end

	return pars_estimate, chi_values

	println("All done!");
end
