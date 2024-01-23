% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;

%% Load data
% Folder containing all one-step datasets
folder = '../Data/OneStep_datasets/';
% One-step metadata file
% Dataset organized by reported latent period in ascending order
master_file = '../Data/OneStep_data.csv';
master_table = readtable(master_file,'ReadVariableNames',true);
files = table2array(master_table(:,'filename'));
xaxis_units = table2array(master_table(:,'xaxis'));
yaxis_units = table2array(master_table(:,'yaxis'));
reportedLP = table2array(master_table(:,'reportedLP'));
host = table2array(master_table(:,'host'));
virus = table2array(master_table(:,'virus'));

%% Plot data
lafigura = figure('Position',[10,10,950,900]);
% Loop through files
for i = 1:numel(files)
    % read file
    file = csvread(strjoin([folder,files(i)],''));
    % plot 
    nexttile;
    scatter(round(file(:,1)),file(:,2),'MarkerFaceColor','k','MarkerEdgeColor','k'); hold on;
    %draw vertical line for reported latent period
    xline(reportedLP(i),'LineStyle',':','LineWidth',2);
    xlabel(xaxis_units(i),'FontSize',18);
    ylabel(yaxis_units(i),'FontSize',18);
    %set x limits and xticks
    xlim([0,ceil(max(file(:,1)))]);
    if max(file(:,1))>=150 %increments of 50
        xticks([0,50:50:ceil(max(file(:,1)))]);
    elseif max(file(:,1))>=60 %increments of 20
        xticks([0,20:20:ceil(max(file(:,1)))]);
    elseif max(file(:,1))>=30 %increments of 10
        xticks([0,10:10:ceil(max(file(:,1)))]);
    else %increments of 2
        xticks([0,2:2:ceil(max(file(:,1)))]);
    end
    xtickangle(0); 
    %set ylims and yticks, 1 per log10
    minlog10 = log10(min(file(:,2)));
    maxlog10 = log10(max(file(:,2)));
    if yaxis_units(i) == "PFU/ml"    
        ylim([10^floor(minlog10),5*10^round(maxlog10)]);
        yticks(10.^(floor(minlog10):round(maxlog10)));
    else %PFU/infected cell
        ylim([min(1,min(file(:,2))),5*10^round(maxlog10)]);
        yticks(10.^(0:round(maxlog10)))
    end
    title([strjoin(['\textit{',host(i),'}'],''),virus(i)],'Interpreter','Latex')

    %aesthetics
    box off;
    set(gca,'YScale','log');
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.03,0.03]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=15;
end

% Add legend to first plot
nexttile(1);
qw{1} = plot(nan, ':','LineWidth',2,'Color','k');
legend([qw{:}], {'First burst'}, 'location', 'southeast');
%%
saveas(gcf,'../Figures/Figure1.svg');