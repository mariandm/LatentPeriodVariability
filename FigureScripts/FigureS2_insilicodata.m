% Code by Marian Dominguez-Mirazo, 2023
% This script plots data previously created for the paraemter inference...
% section of the paper
clear all; close all; clc;

%% Define data ids and simulation id
% List of dataset IDs, see Table 2 main text for parameter details
% See the Parameter Inference folder for details on data generation
% data1: E.coli and lambda
% data2: P.marinus and PHM2
% data3: E.hux and EhV
dataids = ['data1';'data2';'data3'];
% Simulation ids 
ids = 1:10;
%% Other plot information
% Color map
cmap = parula(10);
% Define axis limits for plot
% Xaxis (time, hr) limit for viral dynamics
maxxsV=[12,30,30];
% Xaxis (time, hr) limit to host dynamics (to fill from data)
maxxsH =[];
% Yaxis (virion density, PFU/ml) limit for viral dynamics
maxysV = [5e11,5e10,1e9];
% Yaxis (host density, CFU/ml) limit for host dynamics
maxysH = [1e9,1e9,1e6];
% Lysis rate (hr^-1), to fill from data
etas=[];
%%
% Figure size and position
figure('Position',[10,10,1000,780]);
tcl = tiledlayout(3,3,'TileSpacing','Compact');
% Default to latex interpreter
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% Retrieve original dynamics and plot
% Loop through data ids
for j =1:size(dataids,1)
    dataid = dataids(j,:);
    % Loop through simulation ids
    for i = 1:numel(ids)     
        id=ids(i);
        % Get 'multi-cycle curve' dynamics
        file = strjoin(['../ParameterInference/CreateData/',dataid,'/virOG_',string(id),'.csv'],'');
        vir = readtable(file, 'ReadVariableNames', false);
        vir = table2array(vir);
        % Get original parameters with which data was created
        file = strjoin(['../ParameterInference/CreateData/',dataid,'/parsOG_',string(id),'.csv'],'');
        pars = readtable(file, 'ReadVariableNames', false);
        pars = table2array(pars);
        % Save current lysis rate
        eta = pars(5); % hr^-1
        % Save current coefficient of variation
        cv = pars(6);
        % Define latent period distribution in terms of shape and scale
        shape = 1/(cv^2);
        scale = (cv^2)/eta;
        % Plot latent period distribution
        nexttile(3*(j-1)+1);
        plot(0:0.01:(1/eta)*2,gampdf(0:0.01:(1/eta)*2,shape,scale),...
        'Color',cmap(i,:),'LineWidth',2); hold on;
        % Plot viral dynamics
        nexttile(3*(j-1)+2);
        semilogy(vir(:,1),vir(:,end),'LineWidth',2,...
            'Color',cmap(i,:)); hold on; 
        % Plot total host dynamics
        nexttile(3*(j-1)+3);
        semilogy(vir(:,1),sum(vir(:,2:end-1),2),...
            'LineWidth',2,'Color',cmap(i,:)); hold on; 

    end
    % Store eta per dataset
    etas=[etas,eta];
    % Store xaxis limit per dataset
    maxxsH=[maxxsH,vir(end,1)];
end

%% Aesthetics
% Draw legend using CV values corresponding to simulation ids
% See README in parameter inference for a table
hL = legend('0.5','0.45','0.4','0.35','0.3','0.25','0.2',...
    '0.15','0.1','0.05', 'interpreter','latex');
% Legend title
title(hL,'CV');
% Legend position
hL.Layout.Tile = 'East';
% subplot identifier
texts=['A';'B';'C'];

% Add additional plot features
for j =1:size(dataids,1)
    
    % Distribution plot
    nexttile(3*(j-1)+1);
    % Plot line for latent period mean
    xline(1/etas(j),'--','LineWidth',2.5,'Color','#808080');
    % Aesthetics
    ax=gca;
    ax.FontSize=15;
    box off;
    set(gca,'FontName','Latin Modern Roman')
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    xlabel('Time, hr','FontSize',18);
    ylabel('PDF','FontSize',18);
    % PLot limits
    xlim([0,(1/etas(j))*2]);
    a = xlim;
    b=ylim;
    % Subplot identifier
    text(-0.25*a(2),2.1/2*b(2), texts(j),'FontSize',25,...
        'FontName','Latin Modern Roman');
    % Define xaxis ticks
    xticks([0,0.5/etas(j),1/etas(j),1.5/etas(j),2/etas(j)])
    
    % Viral dynamics plot
    nexttile(3*(j-1)+2);
    % Set yticks, one per log10
    this_ylim = ylim;
    yticks(10.^(floor(log10(this_ylim(1))):floor(log10(this_ylim(2)))));
    % Set y limit
    ylim([10^max([0,floor(log10(this_ylim(1)))]),maxysV(j)]);
    % Set x limit
    xlim([0,maxxsV(j)]);
    %Aesthetics
    ax=gca;
    ax.FontSize=15;
    box off;
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    xlabel('Time, hr','FontSize',18);
    ylabel('Free virus, PFU/ml','FontSize',18);
    
    % Total host plot
    nexttile(3*(j-1)+3);
    % Set yticks, one per log10
    this_ylim = ylim;
    yticks(10.^(0:floor(log10(this_ylim(2)))));
    % Set ylimit
    ylim([1,maxysH(j)]);
    % Set xlimit
    xlim([0,maxxsH(j)]);
    % Aesthetics
    ax=gca;
    ax.FontSize=15;
    box off;
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    xlabel('Time, hr','FontSize',18);
    ylabel('Host, CFU/ml','FontSize',18);
end
%% Save figure
saveas(gcf,'../Figures/FigureS2.svg');