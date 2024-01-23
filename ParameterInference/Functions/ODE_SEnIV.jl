# Code by Marian Dominguez-Mirazo and David Demory, 2023
# See Methods for the mathematical description of the ODE model

function ODE_SEnIV(du,u,p,t)
    mu, K, phi, beta, eta, n = p[1:6]
    n = Int64(n)
    
    # Special case when there are no E compartments
    if n ==0
        S,I,V = u
        du[1] = mu*S*(1-(S+I)/K) - phi*S*V
        du[2] = phi*S*V - eta*I
        du[3] = beta*eta*I - phi*S*V
        
    else
        rate_in = zeros(eltype(phi),n,n)
        S = u[1]
        E = u[2:(n+1)]
        I = u[n+2]
        V = u[n+3]
        
        for i = 1:n
            rate_in[i,i] = (n+1)*eta
        end
        rate_in[1,1] = phi*S*V
        E_in  = [1;E[1:(n-1)]]
        
        du[1] = mu*S*(1-(S+I+sum(E))/K) - phi*S*V
        du[2:(n+1)] = rate_in*E_in - (n+1)*eta*E
        du[n+2] = (n+1)*eta*E[end] - (n+1)*eta*I
        du[n+3] = beta*(n+1)*eta*I - phi*S*V
        
    end
    
    return du
end