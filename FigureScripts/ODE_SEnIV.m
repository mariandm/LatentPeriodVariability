% Code by Marian Dominguez-Mirazo and David Demory, 2023 
% See Methods for the mathematical description of the ODE model


function dy = ODE_SEnIV(t,y,params)

dy = zeros(params.n+3,1);

% Special case when there are no E compartments
if params.n==0
    
    % n = 0
    S = y(1);         % Susceptible cell
    I = y(2);         % Infected
    V = y(3);         % Free virus
    
    %% ODE system
    dy(1) = params.mu*S.*(1-(S+I)./params.K) - params.phi*S.*V;
    dy(2) = params.phi*S.*V - params.eta*I;
    dy(3) = params.beta*params.eta*I - params.phi*S*V;
    
else
    
    % n > 0
    S = y(1);           % Susceptible cell
    E = y(2:params.n+1);       % n Exposed cell classes
    I = y(params.n+2);         % Infected
    V = y(params.n+3);         % Free virus
    
    % funtions transfer
    rate_in = eye(params.n).*repmat((params.n+1)*params.eta,params.n,1);
    rate_in(1,1) = params.phi*S.*V;
    
    E_in  = [1;E(1:params.n-1)];
    
    %% ODE system
    %
    % Susceptible cells
    dy(1) = params.mu*S.*(1-(S+I+sum(E))./params.K) - params.phi*S.*V;
    dy(2:params.n+1) = rate_in*E_in  - (params.n+1)*params.eta*E;
    dy(params.n+2) = (params.n+1)*params.eta*E(end) - (params.n+1)*params.eta*I;
    dy(params.n+3) = params.beta*(params.n+1)*params.eta*I - params.phi*S*V;
end

end
