# Code by Marian  Dominguez-Mirazo, 2023
# This function takes host growth dynamics in the absence of virus and segments the data. 
# The slope is used as a point estimate of the growth rate (mu) 
# The y axis intersection of the plateau section is used as a point estimate for the carrying capacity (K)

function segment_estimate(tsolH,ysolH)
	# Set storage vector
	pars_estimate=zeros(2);
	
	# Estimate growth rate
	# Determine the slope by finding the data section with the highest rate of change
	# Log of the host dynamics
	H = log.(ysolH);
	# Get the difference between time points
	Hdiff = diff(H);
	# Get the floored log10 of the biggest diference in the data
	order = floor(log10.(maximum(Hdiff)));
	# Get all the difference values that are larger than the order of magnitude previously established
	id = findall(x->x>10^order,Hdiff); 
	# Calculate the slope used in that subsection of the data
	data = DataFrame(X=tsolH[id], Y=log.(ysolH[id]));
	mdl = lm(@formula(Y ~ X), data);
	pars_estimate[1] = coef(mdl)[2];

	#Estimate carrying capacity
	pars_estimate[2]  = maximum(ysolH);

	return(pars_estimate)
end