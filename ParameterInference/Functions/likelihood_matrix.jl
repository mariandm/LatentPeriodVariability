# Code by Marian Dominguez-Mirazo, 2023
# This function takes parameter ranges (or assumes a very broad range if none provided) for multiple parameters and calculates a likelihood map for multiple parameter value combinations. 
# It outputs the maximum likelihood parameter combination, a restricted parameter range, and the paraemeter combinations with their corresponding likelihood value.

# Load dependencies
include("../Functions/run_vir.jl")
include("../Functions/get_likelihood.jl")

function likelihood_matrix(params,t_interval_v,Vsol,t_interval_h,Hsol,noise_VirV,noise_VirH,x0,situation,neis,dt,scale)
    # params: Vector of previously known parameters 
    # t_interval_v: Timepoints for viral data
    # Vsol: Viral data
    # t_interval_h: Timepoints for host data
    # Hsol: Total host data
    # noise_VirV: Standard deviation for viral data
    # noise_VirH: Standard deviation for total host data
    # x0: Initial conditions vector
    # situation: [0,1,2,3] describing what parameters to estimate
	    # 0) all parameters are to be estimated
	    # 1) beta is provided, estimate phi, eta, and cv
	    # 2) phi is provided, estimate beta, eta, and cv
	    # 3) beta and phi are provided, estimate eta and cv
    # neis: [[],[],[],[]] neighborhood, range for each parameter, optional
    # dt: Time step
    # scale: Whether the parameter search should be in log scale or linear

    # Get simulation time
    t_interval = [minimum([t_interval_v;t_interval_h]),maximum([t_interval_v;t_interval_h])]

    # Create vectors of values within parameter ranges
    # Parameter range for eta
    # If range is not provided
    if isempty(neis[1]) 	
	eta_max = 1/(5/60);
	eta_min = 1/40;
	vary_eta = 2 .^ LinRange(log2(eta_min), log2(eta_max),10);
    # If range is provided
    else 			
	eta_min = minimum(neis[1]);
        eta_max = maximum(neis[1]);
        
        if scale=="log"
	        vary_eta = 2 .^ LinRange(log2(eta_min), log2(eta_max),10);
        else
        	vary_eta = collect(LinRange(eta_min, eta_max,10));
        end
    end
    vary_eta = unique(round.(vary_eta,digits=3));

    # Parameter range for CV
    # If range is not provided
    if isempty(neis[2]) 	
    	vary_cv = collect(0.5:-0.05:0.05);
    # If range is provided
    else 			
    	cv_min = minimum(neis[2]);
    	cv_max = maximum(neis[2]);
    	vary_cv = collect(LinRange(cv_min,cv_max,40));
    	vary_cv = unique(round.(vary_cv,digits=2));
    end

    # Parameter range for phi
    # If range is not provided
    if isempty(neis[3])
        phi_min = 1e-12;
        phi_max = 1e-6;
        vary_phi = 10 .^ LinRange(log10(phi_min), log10(phi_max),8)
    # If range is provided
    else
        phi_min = minimum(neis[3]);
        phi_max = maximum(neis[3]);
			
	if scale=="log"
       	    vary_phi = 10 .^ LinRange(log10(phi_min),log10(phi_max),10);
	else
	    vary_phi = collect(LinRange(phi_min,phi_max,3)); #10
	end 
    end

    # Parameter range for beta
    # If range is not provided
    if isempty(neis[4])
        beta_min = 1;
        beta_max = 1000;
        vary_beta = 2 .^ LinRange(log2(beta_min), log2(beta_max),8)
    # If range is provided
    else
        beta_min = minimum(neis[4]);
        beta_max = maximum(neis[4]);
                	
        if scale=="log"
       	    vary_beta = 2 .^ LinRange(log2(beta_min), log2(beta_max),10);
        else
            vary_beta = collect(LinRange(beta_min,beta_max,10));
        end
    end
    vary_beta = unique(round.(vary_beta,digits=2));


    # Get parameter value combinations based on which parameters are being estimated
    # All parameters are estimated
    if situation==0 
    	# Get combinations for all parameter vectors		
	vary_matrix = combvec([vary_eta,vary_cv,vary_phi,vary_beta]);
    # Assume burst size is known, estimate phi, eta, cv
    elseif situation==1 
    	vary_matrix = combvec([vary_eta,vary_cv,vary_phi]);
    # Assume phi is known, estimate beta, eta, cv
    elseif situation==2
    	vary_matrix = combvec([vary_eta,vary_cv,vary_beta]);
    # Assume phi and beta are known, estimate eta and cv
    else
    	vary_matrix = combvec([vary_eta,vary_cv]);
    end
    # Number of parameter combinations
    nrows = size(vary_matrix)[1];
    println("Number of parameter combinations to test:");
    println(nrows);
    
    # Create arrays to store simulations
    V_solns_vector = zeros(nrows,length(t_interval_v));
    H_solns_vector = zeros(nrows,length(t_interval_h));
        
    # Run a simulation for each parameter combination
    for i = 1:nrows
       	println(i);
        
        # Assign parameter vector based on which parameters are being estimated
        # All parameters are estimated
    	if situation==0 
            p = [params[1],params[2],vary_matrix[i][3],vary_matrix[i][4],vary_matrix[i][1],vary_matrix[i][2]];
        # Assume burst size is known, estimate phi, eta, cv
        elseif situation==1 
       	    p = [params[1],params[2],vary_matrix[i][3],params[4],vary_matrix[i][1],vary_matrix[i][2]];
       	# Assume phi is known, estimate beta, eta, cv
   	elseif situation==2
   	    p = [params[1],params[2],params[3],vary_matrix[i][3],vary_matrix[i][1],vary_matrix[i][2]];
   	# Assume phi and beta are known, estimate eta and cv
    	else
    	    p = [params[1],params[2],params[3],params[4],vary_matrix[i][1],vary_matrix[i][2]];
    	end

        # Run simulation for the specific parameter combination
        nominales = run_vir(p,x0,t_interval[end],dt);
        # Get simulation indexes that match the total host data timepoints
        idx = indexin(t_interval_h,nominales[:,1]);
	idx = something.(idx);
	# Sum all host compartments to get total host dynamics 
        hnominales = sum(nominales[idx,2:end-1],dims=2);
        # Store simulations
        H_solns_vector[i,:] = hnominales;
        # Get simulation indexes that match the viral data timepoints
        idx = indexin(t_interval_v,nominales[:,1]);
	idx = something.(idx);
	# Store simulations
        vnominales = nominales[idx,end];
        V_solns_vector[i,:] = vnominales;

    end

    # Prevent negative densities and infinite log transformations
    H_solns_vector[H_solns_vector.<0] .= 0; H_solns_vector = H_solns_vector .+1;
    V_solns_vector[V_solns_vector.<0] .= 0; V_solns_vector = V_solns_vector .+1;
	
    # Calculate likelihood 
    likelihood_map = get_likelihood(V_solns_vector, Vsol,H_solns_vector,Hsol,noise_VirV,noise_VirH);

    # Restrict parameter range based on likelihood values
    # Sort parameter combinations by likelihood value and get the 10% top parameter combinations
    sortind = sortperm(likelihood_map[:,1],rev=true)[1:Int64(round(nrows*0.1))];
    top_matrix = mapreduce(permutedims,vcat,vary_matrix[sortind])

    # Borderline cases for new eta range
    # Get the id of the maximum likelihood value for eta
    idx = indexin(top_matrix[1,1],vary_eta); idx=idx[1];
    # if the MLE eta is the first one in the original eta range
    if idx==1
	neis_eta = [vary_eta[idx],vary_eta[idx+1]];
    # if the MLE eta is the last one in the original eta range
    elseif idx==length(vary_eta)
	neis_eta = [vary_eta[idx-1],vary_eta[idx]];
    else 
	neis_eta = [vary_eta[idx-1],vary_eta[idx+1]]
    end

    neis_cv = [minimum(top_matrix[:,2]),maximum(top_matrix[:,2])];

    # Store maximum likelihood point estimates for eta and cv
    LL=zeros(4);
    LL[1] = vary_matrix[sortind[1]][1]; 
    LL[2] = vary_matrix[sortind[1]][2]; 
	
    # All parameters are estimated
    if situation==0 
    	# MLE for phi
	LL[3] = vary_matrix[sortind[1]][3]; 
	# MLE for beta
	LL[4] = vary_matrix[sortind[1]][4]; 
	# New parameter range for phi
	neis_phi = [minimum(top_matrix[:,3]),maximum(top_matrix[:,3])];
	# New parameter range for phi
	neis_beta = [minimum(top_matrix[:,4]),maximum(top_matrix[:,4])];
	neis= [neis_eta,neis_cv, neis_phi, neis_beta];
    # Assume burst size is known, estimate phi, eta, cv
    elseif situation==1
    	# MLE for phi
	LL[3] = vary_matrix[sortind[1]][3];
	# New parameter range for phi
	neis_phi = [minimum(top_matrix[:,3]),maximum(top_matrix[:,3])];
	neis= [neis_eta,neis_cv, neis_phi, []];
    # Assume phi is known, estimate beta, eta, cv
    elseif situation==2 
    	# MLE for beta
	LL[4] = vary_matrix[sortind[1]][3]; #assigning beta
	# New parameter range for beta
	neis_beta = [minimum(top_matrix[:,3]),maximum(top_matrix[:,3])];
	neis= [neis_eta,neis_cv, [], neis_beta];
    # Assume phi and beta are known, estimate eta and cv
    else
	neis= [neis_eta,neis_cv, [], []]
    end

    # Save parameter combinations and corresponsing likelihood values
    ll_values = hcat(mapreduce(permutedims, vcat, vary_matrix),likelihood_map);
	
    return LL, neis, ll_values

end
