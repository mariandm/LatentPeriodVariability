# Code by Marian Dominguez-Mirazo and Jeremy Harris, 2023
# This function calculates the likelihood of multiple simulations given a timeseries dataset of viral and total host dynamics, and the standard error obtained from replicates.
# It returns a vector of likelihoods as large as the number of simulations

function get_likelihood(yfreevirus,Vsol,ytotalhost,HVsol,noiseV,noiseH)

	# yfreevirus: Viral simulations
	# Vsol: Viral data
	# ytotalhost: Total Host from simulations
	# HVsol: Total Host data
	# noiseV: Standard error from viral replicates
	# noiseH: Standard error from total host replicates

	# Number of simulations
	rows = size(yfreevirus,1);
	# Number of viral timepoints
	pts_v = size(yfreevirus,2);
	# Number of host timepoints
	pts_h = size(ytotalhost,2);
	# Set up storage
	LogL = zeros(rows,1);

	#Set normal distribution function to feed data, simulations, and experiment noise
	f(data,pred,noise)=pdf(Normal(pred,noise),data)

	# For each simulation
	for i=1:rows
		logLV = 0;
		logLH = 0;

		#Calculate virus likelihood for each timepoint and add
		for j = 1:pts_v
			logLV = logLV + log(f(log.(Vsol[j]),log.(yfreevirus[i,j]),noiseV[j])+1);
		end

		#Calculate host likelihood for each timepoint and add
		for j = 1:pts_h
			logLH = logLH + log(f(log.(HVsol[j]),log.(ytotalhost[i,j]),noiseH[j])+1);
		end

		#Sum virus and host likelihoods
		LogL[i] = logLV .+ logLH;
	end
	
	return LogL
end
