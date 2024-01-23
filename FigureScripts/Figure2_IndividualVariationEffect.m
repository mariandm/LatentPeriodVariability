% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;

%% Set parameters, see main text Table 1
pars.mu     = 0.1;    % Host growth, hr^-1
pars.K      = 1e9;    % Carrying capacity, CFU/ml
pars.phi    = 1e-7;   % Adsorption rate, ml/(CFUxhr)
pars.beta   = 200;    % Burst size
pars.initS  = 1e8;    % Initial host density, CFU/ml
pars.initV  = 1e6;    % Initial viral density, PFU/ml

%% To generate subplots B and C, set eta and loop through n
pars.eta      = 1/4;    % Lysis rate, hr^-1
ns=[24,100,399];        % Number of E compartments

%% Other options
% Numerical integrator options
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
% Color map
cmap = parula(7);

%% Describe LP distributions in terms of shape and scale
means = 1/pars.eta;
variances = (1/pars.eta)^2 .* (1./(ns+1));
cvs = sqrt(variances)./means;
shapes = means^2./variances;
scales = variances./means;

%% Simulate one-step and plot
figure('Position',[10,10,900,580]);
tcl = tiledlayout(2,2);
cnt = 1;
% Simulate one-step growth curves
for i = 1:numel(ns)
    pars.n = ns(i);
    % Plot distribution
    nexttile(1);
    plot(0:0.01:8,gampdf(0:0.01:8,shapes(i),scales(i)), ...
        'Color',cmap(cnt,:),'LineWidth',2.5); hold on;
    % Simulate the adsorption step
    x0 = zeros(pars.n+3,1);
    x0(1) = pars.initS; x0(end) = pars.initV;
    t = 0:0.01:10/60; % 10 minutes coincubation
    [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
    % Dilute 1000-fold
    x0_2 = ysol(end,:)./1000;
    %Simulate one-step tracking
    t_2 = 0:0.001:6; % 6 hrs after dilution
    [tsol_2,ysol_2] = ode45(@ODE_SEnIV,t_2,x0_2,options,pars);
    % Plot one-step curve
    nexttile(2);
    semilogy(tsol_2,ysol_2(:,end),'Color',cmap(cnt,:),'LineWidth',2.5); hold on;
    % Idenitfy first burst
    id=find(diff(ysol_2(:,end))>0.4,1);
    xline(tsol_2(id),':','LineWidth',3,'Color',cmap(cnt,:)); hold on;
    
    cnt = cnt + 1;
end

%% Aesthetics
nexttile(1);
text(-0.15*8,2.3, 'B','FontSize',25,'FontName','Latin Modern Roman');
title('Latent period distribution','FontSize',18);
xline(4,'--','Color','#808080','LineWidth',2.5);
xlim([0,8]);
box off;
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xlabel('Time, hr','FontSize',18);
ylabel('Probability Distribution','FontSize',18);
% Add legend
xline(9,':','Color','#808080','LineWidth',3) ;
legend('CV = 0.2', 'CV = 0.1', 'CV = 0.05',...
    'mean LP = $T$','First burst',...
    'interpreter','latex','Location','northeast');

nexttile(2);
text(-0.16*6,1.5e6, 'C','FontSize',25,'FontName','Latin Modern Roman');
title('One-step growth curve','FontSize',18);
ylim([1e2, 5e5]);
yticks([1e1 1e2 1e3 1e4 1e5 1e6]);
box off;
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xlabel('Time, hr','FontSize',18);
ylabel('Free virus, PFU/ml','FontSize',18);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To generate subplots D and E, loop through eta and n
etas = [1/9,1/5.9,1/5];

%% Describe LP distributions in terms of shape and scale
means = 1./etas;
variances = (1./etas).^2 .* (1./(ns+1));
cvs = sqrt(variances)./means;
shapes = means.^2./variances;
scales = variances./means;

%% Simulate one-step and plot

% Simulate one-step growth curves
cnt = 1;
for i = 1:numel(ns)
    pars.n = ns(i);
    pars.eta = etas(i);
    % Plot distribution
    nexttile(3);
    plot(0:0.01:18,gampdf(0:0.01:18,shapes(i),scales(i)),...
        'Color',cmap(cnt,:),'LineWidth',2.5); hold on;
    xline(means(i),'--','LineWidth',2.5,'Color',cmap(cnt,:)); hold on;
    % Simulate the adsorption step
    x0 = zeros(pars.n+3,1);
    x0(1) = pars.initS; x0(end) = pars.initV;
    t = 0:0.01:10/60; % 10 minutes
    [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
    % Dilute 1000-fold
    x0_2 = ysol(end,:)./1000;
    t_2 = 0:0.001:9.5;
    [tsol_2,ysol_2] = ode45(@ODE_SEnIV,t_2,x0_2,options,pars);
    % Plot one-step curve
    nexttile(4);
    semilogy(tsol_2,ysol_2(:,end),'Color',cmap(cnt,:),'LineWidth',2.5); 
    hold on;
    % Plot distribution mean
    xline(means(i),'--','LineWidth',2,'Color',cmap(cnt,:)); hold on;
    
    cnt = cnt + 1;
end
% Plot first burst
xline(4,':','LineWidth',3,'Color','#808080'); hold on;

%% Aesthetics
nexttile(3);
text(-0.15*18,1.8, 'D','FontSize',25,'FontName','Latin Modern Roman');
xlim([0,18]);
box off;
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xlabel('Time, hr','FontSize',18);
ylabel('Probability Distribution','FontSize',18);

nexttile(4);
text(-0.16*9.3,1.5e6, 'E','FontSize',25,'FontName','Latin Modern Roman');
xlim([0,9.3]);
ylim([1e2,5e5]);
yticks([1e1 1e2 1e3 1e4 1e5 1e6]);
box off;
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xlabel('Time, hr','FontSize',18);
ylabel('Free virus, PFU/ml','FontSize',18);

%%
saveas(gcf,'../Figures/Figure2.svg');