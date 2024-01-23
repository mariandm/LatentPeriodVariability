% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;

%% Define data ids and simulation ids
% List of dataset IDs, see Table 2 main text for parameter details
% See the Parameter Inference folder for details on data generation
% data1: E.coli and lambda
% data2: P.marinus and PHM2
% data3: E.hux and EhV
dataids = ["data1";"data2";"data3"];
% CV values
cvs = [0.6,0.5:-0.05:0.12];
% Numerical integrator options
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
for i = 1:numel(dataids)
    dataid = dataids(i);
    % Data path
    path = strjoin(['../ParameterInference/CreateData/',dataid,'/'],'');
    first_bursts = [];
    real_cvs = [];
    for j =1:numel(cvs)
        cv = cvs(j);
        % Read parameter values
        file = strjoin([path,'parsOG_1.csv'],'');
        pars_table = readtable(file, 'ReadVariableNames', false);
        pars_table = table2array(pars_table);
        % Read initial conditions
        file = strjoin([path,'x0OG_1.csv'],'');
        init_table = readtable(file, 'ReadVariableNames', false);
        init_table = table2array(init_table);
        % Assign parameters
        pars.mu = pars_table(1);
        pars.K = pars_table(2);
        pars.phi = pars_table(3);
        pars.beta = pars_table(4);
        pars.eta = pars_table(5);
        % Approximate CV (see methods)
        pars.n = round(1./(cv.^2)-1); 
        % Transform back to CV value
        real_cvs = [real_cvs,1 / sqrt(pars.n+1),];
        % Assign initial conditions
        x0 = zeros(pars.n+3,1);
        x0(1) = pars.K/100; x0(end) = pars.K/1000;
        % Simulate coincubation
        t = 0:0.01:10/60; % 10 minutes coincubation
        [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
        % Dilute 1000-fold
        x0_2 = ysol(end,:)./1000;
        %Simulate one-step tracking
        t_2 = 0:0.001:2/pars.eta; % 2 hrs after dilution
        [tsol_2,ysol_2] = ode45(@ODE_SEnIV,t_2,x0_2,options,pars);
        % Get indexes with visible virus
        idx = (ysol_2(:,end) - ysol_2(1,end)) > 0.05;
        % Get corresponding timepoints
        burst_timepoints = tsol_2(idx);
        % Time of first burst
        this_first_burst = burst_timepoints(1);
        first_bursts = [first_bursts, this_first_burst];
    end
    writematrix([real_cvs',first_bursts'],...
        strjoin(['../Data/SimulatedOneSteps/',dataid,'_oneStepTimes.csv'],''),...
        'Delimiter',',')
end