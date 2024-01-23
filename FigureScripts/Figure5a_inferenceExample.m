% Code by Marian Dominguez-Mirazo, 2023
% The script GenerateCIforFig5a.m has been run in advanced. 
% Its output is required to run the present script.
clear all; close all; clc;
%% Define dataset ID and simulation ID
datapath = '../ParameterInference/CreateData';
% Dataset ID
dataid = 'data1';
% Simulation ID
id='6';
% Underlying CV value
cvOG = 0.25;
% Underyingg latent period mean 
TOG = 1;
% Initial conditions, see main text table 2
pars.initS = 1e5;
pars.initV = 1e3;

%% Retrieve simulation replicates

% Retrieve simulated viral replicates
file = strjoin([datapath,'/',dataid,'/vir_noisev_',string(id),'.csv'],'');
virnoise = readtable(file, 'ReadVariableNames', false);
virnoise = table2array(virnoise);
% Calculate replicates mean
mean_virnoise = mean([virnoise(:,2),virnoise(:,3),virnoise(:,4)],2);
% Calculate standard deviation
sd_virnoise = std([virnoise(:,2),virnoise(:,3),virnoise(:,4)],[],2);

% Retrieve simulated total host replicates
file = strjoin([datapath,'/',dataid,'/vir_noiseh_',string(id),'.csv'],'');
hostnoise = readtable(file, 'ReadVariableNames', false);
hostnoise = table2array(hostnoise);
% Calculate replicates mean
mean_hostnoise = mean([hostnoise(:,2),hostnoise(:,3),hostnoise(:,4)],2);
% Calculate standard deviation
sd_hostnoise = std([hostnoise(:,2),hostnoise(:,3),hostnoise(:,4)],[],2);

%% Load host-only MCMC chains
% Read chains
file = strjoin(['../ParameterInference/HostOnlyParams/',dataid,'/hostchain_',string(id),'.csv'],'');
tab = readtable(file, 'ReadVariableNames', false);
tab = table2array(tab);
% Save growth rate chain
mu_chain = tab(:,1);
% Save carrying capacity chain
K_chain = tab(:,2);
% Get point estimates
point_mu = mean(mu_chain);
point_K = mean(K_chain);

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
% Get point estimates
point_ad = mean(ad_chain);
point_beta = mean(beta_chain);
point_eta = mean(eta_chain);
point_cv = mean(cv_chain);

%% Load CI fit
% See script: GenerateCIforFig5a
load(['../Data/vir_CI_',dataid,'_',id,'.mat']);
load(['../Data/host_CI_',dataid,'_',id,'.mat']);

%% Run fit using point estimates
% Numerical integrator options
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
% Run timepoints (hr^-1)
t = 0:0.1:max([virnoise(end,1),hostnoise(end,1)]);
% Set parameters to point estimates
pars.mu = point_mu;
pars.K = point_K;
pars.phi = point_ad;
pars.beta = point_beta;
pars.eta = point_eta;
cv = point_cv;
% Approximate CV (see methods)
pars.n = round(1./(cv.^2)-1);
% Initial conditions
x0 = zeros(pars.n+3,1);
x0(1) = pars.initS; x0(end) = pars.initV;
% Run
[tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
% Save viral dynamics
mean_vir = ysol(:,end);
% Calculate total host dynamics
mean_host = sum(ysol(:,1:end-1),2);

%% Plot 
figure('Position',[10,10,1000,260]);
tcl = tiledlayout(1,3,'TileSpacing','compact');

%% Plot viral dynamics
% Get ranges of CI shaded area
vmin = min(save_vir,[],2);
vmax = max(save_vir,[],2);
nexttile(1);
% Plot CI shaded area
fill([t';flipud(t')],[vmin;flipud(vmax)],'k','FaceColor','#808080',...
    'EdgeColor','none','FaceAlpha',0.5); hold on;
% Plot simulated replicates
errorbar(virnoise(:,1),mean_virnoise,sd_virnoise,'.',...
    'Color','k','LineWidth',1.5);

%% Plot total host dynamics
% Get ranges of CI shaded area
hmin = min(save_host,[],2);
hmin(hmin<0)=0.00001;
hmax = max(save_host,[],2);
hmax(hmax<0)=0.00001;
nexttile(2);
% Plot CI shaded area
fill([t';flipud(t')],[hmin;flipud(hmax)],'k','FaceColor','#808080','EdgeColor','none','FaceAlpha',0.5); hold on;
% Plot simulated replicates
errorbar(hostnoise(:,1),mean_hostnoise,sd_hostnoise,'.','Color','k','LineWidth',1.5);

%% Plot infered distribution
% Get CV 95% CI chain steps
cv_95chain = cv_chain(cv_chain > quantile(cv_chain,0.025) ...
    & cv_chain < quantile(cv_chain,0.975));
% Get lysis rate 95% CI chain steps
eta_95chain = eta_chain(eta_chain > quantile(eta_chain,0.025) ...
    & eta_chain < quantile(eta_chain,0.975));
% Shuffle CV chain
cv_95_shuffle = cv_95chain(randperm(length(cv_95chain)));
% Shuffle lysis rate chain
eta_95_shuffle = eta_95chain(randperm(length(eta_95chain)));
% Decribe LP distributions in terms of shape and scale
shapes = 1 ./ (cv_95_shuffle .^2);
scales = (cv_95_shuffle .^2) ./  eta_95_shuffle;
% Get chain latent period means
Ts = 1 ./ eta_95_shuffle;
% Latent period mean range
maxT = max(Ts);
minT = min(Ts);
% Set x axis for latent period distribution
xs = 0:0.01:maxT*2;
% Get latent period distributions
gammas = zeros(numel(Ts),numel(xs));
for j = 1:numel(Ts)
    gammas(j,:) = gampdf(xs,shapes(j),scales(j));
end
% Get latent period distribution ranges across x axis
maxg = max(gammas,[],1);
ming = min(gammas,[],1);
% Calculate underlying latent period distribution
shape = 1/cvOG^2;
scale = cvOG^2 / TOG;
mean_shape = 1/point_cv^2;
mean_scale = point_cv^2 / point_eta;

% Plot latent period distribution
% This is ugly, please close your eyes
nexttile(3);
% First, plotting in the order we want the legend to appear:
% Plot underlying distribution
plot(0:0.01:maxT*2,gampdf(0:0.01:maxT*2,shape,scale),...
    'Color','k','LineWidth',1.5); hold on;
% Plot underlying mean(LP)
xline(TOG,'--','Color','k','LineWidth',1.5);
% Plot 95% CI as a geometric figure
% Confidence interval for mean latent period
lelimit = ylim;
fill([min(1./eta_95chain),min(1./eta_95chain),...
    max(1./eta_95chain),max(1./eta_95chain)],...
[lelimit(1),lelimit(2),lelimit(2),lelimit(1)],'k',...
'FaceColor','#808080','EdgeColor','#808080','FaceAlpha',0.5); hold on;
% Confidence interval for gamma distribution
fill([xs';flipud(xs')],[ming';flipud(maxg')],'k',...
    'FaceColor','#808080','EdgeColor','#808080','FaceAlpha',0.5); hold on;

%Plot in the order we want the lines to appear (above the fill):
%Plot underlying distribution
plot(0:0.01:maxT*2,gampdf(0:0.01:maxT*2,shape,scale),...
    'Color','k','LineWidth',1.5);
%PLot underlying mean(LP)
xline(TOG,'--','Color','k','LineWidth',1.5);

%% ~aesthetics~

% Free virus dynamics
nexttile(1);
text(-12,3e8, 'A','FontSize',35,'FontName','Latin Modern Roman');
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
set(gca, 'YScale', 'log');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
yticks([1e3 1e4 1e5 1e6 1e7 1e8 1e9]);
box off;
xlabel('Time, hr','FontSize',18);
ylabel('PFU/ml','FontSize',18);
title('Free virus','FontSize',18);
xlim([0,12]);

% Host dynamics
nexttile(2);
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman')
set(gca, 'YScale', 'log');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
yticks([1e2 1e3 1e4 1e5 1e6 1e7 1e8]);
box off;
ylim([10,1e9]);
xlim([0,12]);
xlabel('Time, hr','FontSize',18);
ylabel('CFU/ml','FontSize',18);
title('Host','FontSize',18);

% Ditribution
nexttile(3);
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
yticks(0:0.5:2);
box off;
xlabel('Time, hr','FontSize',18);
ylabel('Probability Distribution','FontSize',18);
title('Probability Distribution','FontSize',18);
hL = legend('Underlying distribution', 'Underlying mean (LP)','95\% CI prediction','interpreter','latex');
hL.Layout.Tile = 'East';
xlim([0,maxT*2]);

%% Save figure
saveas(gcf,'../Figures/Figure5a.png');
