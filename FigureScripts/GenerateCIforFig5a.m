% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;

%%
datapath = '../ParameterInference/CreateData/';
dataid = 'data1';
id='6';
pars.initS = 1e5;
pars.initV = 1e3;
%% Load simulated replicates (just to get experiment time)
% Viral replicates
file = strjoin([datapath,'/',dataid,'/vir_noisev_',string(id),'.csv'],'');
virnoise = readtable(file, 'ReadVariableNames', false);
virnoise = table2array(virnoise);
% Total host replicates
file = strjoin([datapath,'/',dataid,'/vir_noiseh_',string(id),'.csv'],'');
hostnoise = readtable(file, 'ReadVariableNames', false);
hostnoise = table2array(hostnoise);

%% Load host-only MCMC chains
% Read chains
file = strjoin(['../ParameterInference/HostOnlyParams/',dataid,'/hostchain_',string(id),'.csv'],'');
tab = readtable(file, 'ReadVariableNames', false);
tab = table2array(tab);
% Save growth rate chain
mu_chain = tab(:,1);
% Save carrying capacity chain
K_chain = tab(:,2);

%% Load virus-host MCMC chains
% Read chains
file = strjoin(['../ParameterInference/VirusHostParams/step2_MCMC/round2/',dataid,'/viruschain_',string(id),'_round2.csv'],'');
tab = readtable(file, 'ReadVariableNames', false);
tab = table2array(tab);
% Save adsorption rate chain
ad_chain = 10.^tab(:,1);
% Save burst size chain
beta_chain = tab(:,2);
% Save lysis rate chain
eta_chain = tab(:,3);
% Save CV chain
cv_chain = tab(:,4);


%% Do 1000 chain value combinations
nruns=1000;
% Sample indexes
s1 = randsample(1:numel(mu_chain),nruns);
s2 = randsample(1:numel(K_chain),nruns);
s3 = randsample(1:numel(ad_chain),nruns);
s4 = randsample(1:numel(beta_chain),nruns);
s5 = randsample(1:numel(eta_chain),nruns);
s6 = randsample(1:numel(cv_chain),nruns);
% Numerical integration options
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
% Simulation time
t = 0:0.1:max([virnoise(end,1),hostnoise(end,1)]);

% Run nruns simulations with shuffled combinations of parameters
for i=1:nruns
    pars.mu = mu_chain(s1(i));
    pars.K = K_chain(s1(i));
    pars.phi = ad_chain(s3(i));
    pars.beta = beta_chain(s4(i));
    pars.eta = eta_chain(s5(i));
    cv = cv_chain(s6(i));
    % Approximate CV (see methods)
    pars.n = round(1./(cv.^2)-1);
    % Initial conditions
    x0 = zeros(pars.n+3,1);
    x0(1) = pars.initS; x0(end) = pars.initV;
    % Run 
    [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
    % Save viral dynamics
    save_vir(:,i) = ysol(:,end);
    % Calculate total host dynamics
    save_host(:,i) = sum(ysol(:,1:end-1),2);
end
%% Save simulations
save(['../Data/vir_CI_',dataid,'_',id,'.mat'],'save_vir');
save(['../Data/host_CI_',dataid,'_',id,'.mat'],'save_host');