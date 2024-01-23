% Code by Marian Dominguez-Mirazo, 2023
% The script GenerateOneStepTimesforFig5b.m has been run in advanced. 
% Its output is required to run the present script.
clear all; close all; clc;

%% Define data ids and simulation id correspondence to CV
% List of dataset IDs, see Table 2 main text for parameter details
% See the Parameter Inference folder for details on data generation
% data1: E.coli and lambda
% data2: P.marinus and PHM2
% data3: E.hux and EhV
dataids = ["data1";"data2";"data3"];
% Corresponding Latent Period mean
true_Ts = [1,5,6];
% Simulation ids 
ids = [1,2,3,4,5,6,7,8];
% Corresponding CV values
cvs = [0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15];
%%
figure('Position',[10,10,800,800]);
tcl = tiledlayout(3,2,'TileSpacing','compact');
%% Retrieve Bayesian parameter estimates and calculate confidence intervals
% Loop through dataset
for j =1:size(dataids,1)
    % Current dataset ID
    dataid=dataids(j,:);
    % Current Latent period mean
    true_T = true_Ts(j);
    % Create confidence interval storage
    errorbar_eta = zeros(numel(ids),3);
    errorbar_cv = zeros(numel(ids),3);
    % Loop through simulation ID
    for i = 1:numel(ids)
        % Current simulation ID
        id = ids(i);
        % Corresponding CV
        cv = cvs(i);
        % Get MCMCM chains
        file = strjoin(['../ParameterInference/VirusHostParams/step2_MCMC/round2/',dataid,'/viruschain_',string(id),'_round2.csv'],'');
        tab = readtable(file, 'ReadVariableNames', false);
        tab = table2array(tab);
        % Lysis rate chain
        eta_chain = tab(:,3);
        % CV chain
        cv_chain = tab(:,4);
        
        %% Calculate Confidence Intervals
        % Point estimate
        point_eta = mean(eta_chain);
        point_cv = mean(cv_chain);
        % 95% CI values for lysis rate
        q25_eta = quantile(eta_chain,0.025);
        q975_eta = quantile(eta_chain,0.975);
        % 95% CI values for CV
        q25_cv = quantile(cv_chain,0.025);
        q975_cv = quantile(cv_chain,0.975);
        % Save info in array
        % Save latent period mean point estimate
        errorbar_eta(i,1) = 1/point_eta;
        % Save Latent period mean error bar (from confidence interval)
        errorbar_eta(i,2:3) = [abs(1/q975_eta-1/point_eta),abs(1/q25_eta-1/point_eta)];
        % Save CV point estimate
        errorbar_cv(i,1) = point_cv;
        % Save CV error bar (from confidence interval)
        errorbar_cv(i,2:3) = [abs(q25_cv-point_cv),abs(q975_cv-point_cv)];
    end
    
    %% Read first burst times from simulated one-step data
    file = strjoin(['../Data/SimulatedOneSteps/',dataid,'_oneStepTimes.csv'],'');
    onestep = readtable(file, 'ReadVariableNames', false);
    onestep = table2array(onestep);

    %% Plot 
    
    % Latent period mean parameter inference
    nexttile(2*(j-1)+1);
    % Plot the underlying latent period mean (with which data was created)
    plot([0,1.05],[true_T,true_T],'--','LineWidth',1.5,'Color','k'); 
    hold on;
    % Plot the errorbar per simulation
    errorbar(cvs,errorbar_eta(:,1),errorbar_eta(:,2),errorbar_eta(:,3),...
        '.','Color','#707070','LineWidth',2.5); hold on;
    % Plot the point estimates
    scatter(cvs,errorbar_eta(:,1),'r','filled','SizeData',50);
    % PLot the first burst size prediction
    plot([onestep(:,1);0],[onestep(:,2);true_T],':','Color','#808080','LineWidth',3);
    % Aesthetics
    ylabel('Estimated T, hr','interpreter','latex','FontSize',18);
    ylim([0,1.5*true_T]);
    xlim([0,0.55]);
    xticks(0:0.1:0.5);
    box off;
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=15;
    pbaspect([1.5 1 1])
    
    % CV parameter inference
    nexttile(2*(j-1)+2);
    % Plot the CV (with which data was created)
    plot([0,1.05],[0,1.05],'--','LineWidth',1.5,'Color','k'); hold on;
    % Plot the errorbar per simulation
    errorbar(cvs,errorbar_cv(:,1),errorbar_cv(:,2),errorbar_cv(:,3),...
        '.','Color','#707070','LineWidth',2.5); hold on;
    % Plot the point estimates
    scatter(cvs,errorbar_cv(:,1),'r','filled','SizeData',50);
    % Aesthetics
    ylabel('Estimated CV','interpreter','latex','FontSize',18);
    xlim([0,0.55]);
    ylim([0,0.55]);
    xticks(0:0.1:0.5);
    yticks(0:0.1:0.5);
    box off;
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=15;
    pbaspect([1 1 1])
end

% Add x axis label to the bottom plots only
nexttile(5);
xlabel('Original CV','interpreter','latex','FontSize',18);
% Include a single legend for all panels
hL = legend('Underlying value','95% CI',...
    'Point Estimate','One-step estimate');
hL.Layout.Tile = 'East';

nexttile(6);
xlabel('Original CV','interpreter','latex','FontSize',18);

%% Save figure
saveas(gcf,'../Figures/Figure5b.png');