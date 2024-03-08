%% Create the 2x1 subplot for Frontiers Papers
clear all; close all; clc;

%% Import Functions
addpath('../functions')

% Load & Save Results of the Surrogates for Significance and Non Linear
% Dynamics

% Significance
load('IS_HH_Sig.mat');
Significant_Surr_Results=sum(knnSE_significance,1);

% Distribution of the IS Values
delta_x = delta;
IS_Values=knnSE;

% Non Linear Dynamics
load('IS_HH_IAFFT.mat');
NonLinear_Surr_Results=sum(knnSE_significance,1);

% Clear all variables execpt the results
clearvars -except Significant_Surr_Results IS_Values NonLinear_Surr_Results delta_x

% Plot
fig=figure;
fig.WindowState='maximized';
subplot_tight(1,2,1,[0.13 0.06]);
errorbar(delta_x, median(IS_Values), median(IS_Values)-quantile(IS_Values,0.05,1),quantile(IS_Values,0.95,1)-median(IS_Values),'LineWidth',2,'Color','k');
ylabel('S_Y','FontWeight','bold');
ylim([-0.1 2.6]);
xlabel('\sigma','FontWeight','bold');
xticks(delta_x);
xlim([-0.1 3.1]);
ax=gca;
ax.FontSize=22;

subplot_tight(1,2,2,[0.13 0.06]);
Data=[Significant_Surr_Results' NonLinear_Surr_Results'];
bar(delta_x,Data);
ylim([0 115]);
xlabel('\sigma','FontSize',25,'FontWeight','bold');
xticks(delta_x);
xlim([-0.1 3.1]);
ylabel('%','FontSize',25,'FontWeight','bold');
yline(5,'-','LineWidth',2,'Color','r');
yticks(0:20:100);
% title('$$linear\ Model: Information\ Storage\ significance - removed\ nonlinearity$$','Interpreter','latex','FontSize',25)
legend({'Presence of Self-Dependencies' 'Presence of Nonlinearities'},'Box','off','Location','north','NumColumns',2);
ax=gca;
ax.FontSize=22;
AddLetters2Plots(fig,{'(a)','(b)'} ,'HShift',[0,0], 'VShift', 0, 'Direction', 'LeftRight','FontSize',20);

% Save Figure
% exportgraphics(fig,'IS_HenonMap_KNN_errorbar_bar_new.eps','BackgroundColor','none');
