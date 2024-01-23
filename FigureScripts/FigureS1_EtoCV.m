% Code by Marian Dominguez-Mirazo, 2023
clear all; close all; clc;

%% Set interpreter 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
%% Plot E to CV relation
ns = 0:99;
cvs = 1./sqrt(ns+1);
figure('Position',[10,10,400,300]);
plot(ns,cvs,'LineWidth',2,'Color','k');

%% Aesthetics
xlabel('Number of E compartments','FontSize',18);
ylabel('Coefficient of Variation','FontSize',18);
xticks([0,25:25:100]);
box off;
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.03,0.03]);
set(gca,'LineWidth',0.6);

%%
saveas(gcf,'../Figures/FigureS1.svg');