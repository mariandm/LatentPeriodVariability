% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;
%% Define ids and underlying parameters
% List of dataset IDs, see Table 2 main text for parameter details
% See the Parameter Inference folder for details on data generation
% data1: E.coli and lambda
% data2: P.marinus and PHM2
% data3: E.hux and EhV
dataids = ["data1";"data2";"data3"];
% Corresponding adsorption rate
true_ads = [1e-8,9.3e-10,1.5e-7];
% Corresponding burst size
true_betas = [200,40,800];
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
    % Current adsorption rate
    true_ad = true_ads(j);
    % Current burst size
    true_beta = true_betas(j);
    % Create confidence interval storage
    errorbar_ad = zeros(numel(ids),3);
    errorbar_beta = zeros(numel(ids),3);
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
        % Adsorption rate chain
        ad_chain = 10.^tab(:,1);
        % Burst size chain
        beta_chain = tab(:,2);
        
        %% Calculate Confidence Intervals
        % Point estimates
        point_ad = mean(ad_chain);
        point_beta = mean(beta_chain);
        % 95% CI values for adsorption rate
        q25_ad = quantile(ad_chain,0.025);
        q975_ad = quantile(ad_chain,0.975);
        % 95% CI values for burst size
        q25_beta = quantile(beta_chain,0.025);
        q975_beta = quantile(beta_chain,0.975);
        %set in adecuate structure for plotting
        % Save info in array
        % Save adsorption rate point estimate
        errorbar_ad(i,1) = point_ad;
        % Save adsorption rate error bar (from confidence interval)
        errorbar_ad(i,2:3) = [abs(q25_ad-point_ad),abs(q975_ad-point_ad)];
        % Save burst size point estimate
        errorbar_beta(i,1) = point_beta;
        % Save burst size error bar (from confidence interval)
        errorbar_beta(i,2:3) = [abs(q25_beta-point_beta),abs(q975_beta-point_beta)];
    end

    %% Plot 
    
    % Adsorption rate parameter inference
    nexttile(2*(j-1)+1);
    % Plot the errorbar per simulation
    errorbar(cvs,errorbar_ad(:,1),errorbar_ad(:,2),errorbar_ad(:,3),...
        '.','Color','k','LineWidth',1.5); hold on;
    % Plot the point estimates
    scatter(cvs,errorbar_ad(:,1),'r','filled','SizeData',50);
    % Plot the underlying adsorption rate (with which data was created)
    yline(true_ad,'--','LineWidth',2.5,'Color','#808080');
    % Aesthetics
    ylabel({'Estimated adsorption rate ($\phi$)';'ml/(CFU$\times$hr)'},'interpreter','latex','FontSize',18);
    ylim([10^(floor(log10(true_ad))-1),10^(ceil(log10(true_ad))+1)]);
    xlim([0,0.55]);
    xticks(0:0.1:0.5);
    box off;
    set(gca, 'YScale', 'log')
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=15;
    pbaspect([1.5 1 1])
    
    % Burst size parameter inference
    nexttile(2*(j-1)+2);
    % Plot the errorbar per simulation
    errorbar(cvs,errorbar_beta(:,1),errorbar_beta(:,2),errorbar_beta(:,3),'.','Color','k','LineWidth',1.5); hold on;
    % Plot the point estimates
    scatter(cvs,errorbar_beta(:,1),'r','filled','SizeData',50);
    % Plot the underlying burst size (with which data was created)
    yline(true_beta,'--','LineWidth',2.5,'Color','#808080');
    % Aesthetics
    ylabel('Estimated Burst Size ($\beta$)','interpreter','latex','FontSize',18);
    ylim([0,2*true_beta]);
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
end

% Add x axis label to the bottom plots only
nexttile(5);
xlabel('Original CV','interpreter','latex','FontSize',18);
% Include a single legend for all panels
hL = legend('95% CI','MLE','Underlying value','One-step estimate');
hL.Layout.Tile = 'East';

nexttile(6);
xlabel('Original CV','interpreter','latex','FontSize',18);

%% Save figure
saveas(gcf,'../Figures/FigureS3.svg');