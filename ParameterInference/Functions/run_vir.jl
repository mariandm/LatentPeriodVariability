# Code by Marian Dominguez-Mirazo, 2023
# This function takes a vector of parameters, initial conditions, final run time and time step, and runs the model approximating the CV in the parameter vector to the closest possible CV value described using an n integer number of E compartments. 

# Load dependencies
include("../Functions/ODE_SEnIV.jl")

function run_vir(p,x0,tf,dt)
    # Find the closest integer n to decribe the parameter CV
    cv = p[6];
    n = 1/cv^2 -1;
    round_n = Int64(round(n));
    params = p[1:5]; push!(params,round_n);

    # Initial conditions
    u0=zeros(round_n+3,1);
    u0[1] = x0[1];
    u0[end] = x0[2];

    # Run simulation
    tspan = (0,tf);
    prob = ODEProblem(ODE_SEnIV,u0,tspan,params);
    y_vir = solve(prob, Tsit5(),saveat=dt,abstol = 1e-8, reltol =1e-8);

    # Make solver output into an array
    ymod_vir = reverse(rotl90(hcat(y_vir.u...)),dims=1);
    # Add time in the first column
    ymod_vir = hcat(y_vir.t,ymod_vir);

    return(ymod_vir)
end
