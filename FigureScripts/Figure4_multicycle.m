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

%% Parameters for reference one-step
reference_mean = 4;
reference_cv = 0.15;
reference_n = round(1/reference_cv^2-1);

%% Parameters for relevant param combinations
colors = ["#4DBEEE","red","#69cf3b"];
la_mean = [3.4,reference_mean,4.6];
la_eta = 1./la_mean;
la_cv = [0.09,reference_cv,0.21];
la_n = round(1./(la_cv.^2)-1);

%% Plot B, C: Multi-cycle response curves
figure('Position',[10,10,850,300]);
tcl = tiledlayout(1,2);

for i=1:numel(la_mean)
    % Set parameters
    pars.eta = la_eta(i);
    pars.n = la_n(i);
    % Initial conditions
    x0 = zeros(pars.n+3,1);
    x0(1) = pars.initS; x0(end) = pars.initV;
    % Simulate nulti-cycle response curve for 18hrs
    t = 0:0.01:18;
    [tsol,ysol] = ode45(@ODE_SEnIV,t,x0,options,pars);
    % Plot viral dynamics
    nexttile(1);
    semilogy(tsol,ysol(:,end),'Color',colors(i),'LineWidth',2.5); 
    hold on;
    % Plot total host dynamics
    nexttile(2);
    semilogy(tsol,sum(ysol(:,1:end-1),2),'Color',colors(i),'LineWidth',2.5); 
    hold on;
    if i==2 %save the reference run
        tsol_reference = tsol;
        vsol_reference = ysol(:,end);
        hsol_reference = sum(ysol(:,1:end-1),2);
    end
end

%%
% Viral dynamics: Represent noise in experiments
nexttile(1);
fill([tsol_reference;flipud(tsol_reference)],...
    [(vsol_reference)*1.4; flipud((vsol_reference)*0.6)],...
    'r','FaceAlpha',0.3)

% Total host dynamics: Represent noise in experiments
nexttile(2);
fill([tsol_reference;flipud(tsol_reference)],...
    [(hsol_reference+0.05)*1.4;flipud(hsol_reference+0.05)*0.6],...
    'r','FaceAlpha',0.3)

% Replot the reference for it to appear in the front
nexttile(1);
semilogy(tsol_reference,vsol_reference,'Color','red','LineWidth',2.5);
nexttile(2);
semilogy(tsol_reference,hsol_reference,'Color','red','LineWidth',2.5);

%% Aesthetics
nexttile(1);
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
xlabel('Time, hr','FontSize',18);
ylabel('PFU/ml','FontSize',18);
xlim([0,18]);
ylim([1 8e10])
yticks([1e0 1e2 1e4 1e6 1e8 1e10]);
title('Free virus');
text(-4,3e11, 'B','FontSize',25,'interpreter','latex');
legend('mean LP = 3.4 hr, CV = 0.09',...
    'mean LP = 4 hr, CV = 0.15',...
    'mean LP = 4.6 hr, CV = 0.21',...
    'Experiment noise', ...
    'Location', 'SouthEast')
box off;
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

nexttile(2);
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman');
xlabel('Time, hr','FontSize',18);
ylabel('CFU/ml','FontSize',18);
yticks([1e0 1e2 1e4 1e6 1e8 1e10]);
xlim([0,18]);
ylim([1 1e10])
title('Total host');
text(-4,6e10, 'C','FontSize',25,'interpreter','latex');
box off;
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

%%
saveas(gcf,'../Figures/Figure4.svg');