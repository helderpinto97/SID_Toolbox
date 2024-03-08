clear all; close all; clc;

%% Import Functions
addpath('functions')

% Load & Save Results of the Surrogates for Significance and Non Linear
% Dynamics

% Significance
load('MIR_VAR_Significance.mat');
tmp_prctile_knn = prctile(knnMIR_surr_varC,95,3);
knnSE_significance = knnMIR_varC > tmp_prctile_knn;
Significant_Surr_Results=sum(knnSE_significance,2);

% Distribution MIR Values - KNN
x = C_arr;
MIR_Values=knnMIR_varC;

% Distribution MIR Theoretical Values
MIR_Values_Th=linMIR_theoretical_varC;

% Non Linear Dynamics
load('MIR_VAR_IAFFT.mat');
tmp_prctile_knn = prctile(knnMIR_surr_varC,95,3);
knnSE_significance = knnMIR_varC > tmp_prctile_knn;
NonLinear_Surr_Results=sum(knnSE_significance,2);

% Clear all variables execpt the results
clearvars -except Significant_Surr_Results MIR_Values NonLinear_Surr_Results x MIR_Values_Th

% Plot
% HTML Gray
htmlGray = [128 128 128]/255;
% Figure
fig=figure;
fig.WindowState='maximized';
subplot_tight(1,2,1,[0.13 0.07]);
errorbar(x, median(MIR_Values,2), median(MIR_Values,2)-quantile(MIR_Values,0.05,2),quantile(MIR_Values,0.95,2)-median(MIR_Values,2),'LineWidth',2,'Color','k','DisplayName','KNN');
hold on;
plot(x,MIR_Values_Th,'LineWidth',2,'Color',htmlGray,'DisplayName','Theoretical Value');
hold off;
legend('Box','off');
ylabel('I_{X;Y}','FontWeight','bold');
% ylim([-0.1 2.6]);
xlabel('$C$','FontWeight','bold','Interpreter','latex');
xticks(x);
xlim([-0.05 1.05]);
ax=gca;
ax.FontSize=21;

subplot_tight(1,2,2,[0.13 0.07]);
Data=[Significant_Surr_Results NonLinear_Surr_Results];
bar(x,Data);
ylim([0 115]);
xlabel('$C$','FontSize',25,'FontWeight','bold','Interpreter','latex');
xticks(x);
xlim([-0.05 1.05]);
ylabel('%','FontSize',25,'FontWeight','bold');
yline(5,'-','LineWidth',2,'Color','r');
yticks(0:20:100);
% title('$$linear\ Model: Information\ Storage\ significance - removed\ nonlinearity$$','Interpreter','latex','FontSize',25)
legend({'Presence of Coupling' 'Presence of Nonlinearities'},'Box','off','Location','north','NumColumns',2);
ax=gca;
ax.FontSize=21;
AddLetters2Plots(fig,{'(a)','(b)'} ,'HShift',[0,0], 'VShift', 0, 'Direction', 'LeftRight','FontSize',20);

% Save Figure
% exportgraphics(fig,'MIR_VAR_KNN_errorbar_bar_new.eps','BackgroundColor','none');

