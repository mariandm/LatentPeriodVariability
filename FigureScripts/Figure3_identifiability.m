% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;

%% Set parameters, see main text Table 1
pars.mu     = 0.1;    % Host growth, hr^-1
pars.K      = 1e9;    % Carrying capacity, CFU/ml
pars.phi    = 1e-7;   % Adsorption rate, ml/(CFUxhr)
pars.beta   = 200;    % Burst size
pars.initS  = 1e8;    % Initial host density, CFU/ml
pars.initV  = 1e6;    % Initial viral density, PFU/ml

% Numerical integrator options
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
% Colormap
cmap = parula(20);

%% To generate Figure 3
% Loop through combinations of eta and n values,
% simulate one-step curves and calculate the error between them.

% CV values
cvs=0.05:0.01:0.25;
% n can only be integer numbers due to model structure, see Methods, FS1
ns = 1./cvs.^2-1;
ns = round(ns); 
% Get approximated cv
cvs = 1./sqrt(ns+1);
% Eta values
etas = 1./(2:0.2:6);

%% Simulate reference one-step

% Set parameters
reference_mean = 4;
reference_cv = 0.15;
pars.eta = 1/reference_mean;
reference_n = round(1/reference_cv^2-1);
pars.n = reference_n;
% Initial conditions
x0 = zeros(pars.n+3,1);
x0(1) = pars.initS; x0(end) = pars.initV;
% Simulate coincubation for 10 minutes
t = 0:0.01:10/60;
[tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
% Dilute 1000-fold
x0_2 = ysol(end,:)./1000;
% Simulate one-step tracking
t_2 = 0:0.001:6;
[tsol_reference,ysol_reference] = ode45(@ODE_SEnIV,t_2,x0_2,options,pars);

%% Simulate one-steps for multiple eta and cv combinations
% Create storage for viral dynamics
logsquares = zeros(numel(ns),numel(etas));
saving_dynamics = [];

for j = 1:numel(etas)
    pars.eta = etas(j);
    for i = 1:numel(ns)
        pars.n = ns(i);
        % Initial conditions
        x0 = zeros(pars.n+3,1);
        x0(1) = pars.initS; x0(end) = pars.initV;
        % Simulate adsorption step for 10 minutes
        t = 0:0.01:10/60;
        [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
        % Dilute 1000-fold
        x0_2 = ysol(end,:)./1000;
        % Simulate one-step tracking
        t_2 = 0:0.001:6;
        [tsol_2,ysol_2] = ode45(@ODE_SEnIV,t_2,x0_2,options,pars);
        %Calculate squared error
        logsquares(i,j) = sum((log10(ysol_reference(:,end)) - ...
            log10(ysol_2(:,end))).^2);
        % Save viral dynamics
        saving_dynamics = [saving_dynamics,ysol_2(:,end)];
    end
end

%% Plot
figure('Position',[10,10,1200,300]);
tcl = tiledlayout(1,3);

%% A: Plot squared error grid
nexttile(1);
imagesc(1./etas,cvs,log10(logsquares)); hold on;
colormap gray;
xlabel('Mean latent period, hr','interpreter','latex','FontSize',18);
ylabel('Coefficient of Variation','interpreter','latex','FontSize',18);
xticks(2:6);
a = colorbar;
set(gca,'Ydir','normal');
set(gca,'FontName','Latin Modern Roman');
ax=gca;
ax.FontSize=15;
title('Log10 squared error','FontSize',18,'interpreter','latex');
pbaspect([1 1 1])
%% A: Add marks to relevant value combinations
colors = ["#4DBEEE","red","#69cf3b"];
la_mean = [3.4,reference_mean,4.6];
la_eta = 1./la_mean;
la_cv = [0.09,reference_cv,0.21];
la_n = round(1./(la_cv.^2)-1);
for i=1:numel(la_mean)
    text(la_mean(i)-0.05,la_cv(i), 'X','Color', colors(i),...
        'FontSize',12,'fontweight','bold');
end
text(0.5,0.27, 'A','FontSize',25,'interpreter','latex');

%% B: Plot one-step growth curves 
nexttile(2);
% Reference first for legend purposes
semilogy(tsol_reference,ysol_reference(:,end),...
    'Color','r','LineWidth',2.5); hold on;
% Add shaded region to represent noise in experiments
fill([tsol_reference;flipud(tsol_reference)],...
    [(ysol_reference(:,end))*1.4;flipud((ysol_reference(:,end))*0.6)],...
    'r','FaceAlpha',0.3)

for i=[1,3]
    %Run examples
    pars.eta = la_eta(i);
    pars.n = la_n(i);
    % Initial conditions
    x0 = zeros(pars.n+3,1);
    x0(1) = pars.initS; x0(end) = pars.initV;
    % Simulate adsorption step for 1 minute
    t = 0:0.01:10/60;
    [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
    % Dilute 1000-fold
    x0_2 = ysol(end,:)./1000;
    % Simulate one-step tracking
    t_2 = 0:0.001:6;
    [tsol_2,ysol_2] = ode45(@ODE_SEnIV,t_2,x0_2,options,pars);
    %Plot examples
    semilogy(tsol_2,ysol_2(:,end),'Color',colors(i),'LineWidth',2.5);
    hold on;
end

%Plot reference (again to have it on front)
semilogy(tsol_reference,ysol_reference(:,end),'Color','r','LineWidth',2.5);

%Aesthetics
yticks([1e1 1e2 1e3 1e4 1e5 1e6]);
ax=gca;
ax.FontSize=15;
legend('Reference one-step','Experiment noise','Location','NorthWest');
set(gca,'FontName','Latin Modern Roman');
xlabel('Time, hr','FontSize',18);
ylabel('Free virus, PFU/ml','FontSize',18);
title('One-step growth curve');
text(-1.5,2e6, 'B','FontSize',25,'interpreter','latex');
box off;
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

%% C: Plot latent period distributions
nexttile(3);
% Write distributions in terms of shape and scale
means = 1./la_eta;
variances = (1./la_eta).^2 .* (1./(la_n+1));
shapes = means.^2./variances;
scales = variances./means;
for i =1:numel(la_eta)
    plot(0:0.01:8,gampdf(0:0.01:8,shapes(i),scales(i)),...
        'Color',colors(i),'LineWidth',2.5); 
    hold on;
end

% Aesthetics
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
xlabel('Time, hr','FontSize',18);
ylabel('Probability distribution','FontSize',18);
title('Latent period distribution');
text(-2,2.15, 'C','FontSize',25,'interpreter','latex');
box off;
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
legend('mean LP = 3.4 hr, CV = 0.09',...
    'mean LP = 4 hr, CV = 0.15',...
    'mean LP = 4.6 hr, CV = 0.21',...
    'Location', 'North')
ylim([0,2]);
%%
saveas(gcf,'../Figures/Figure3.svg');