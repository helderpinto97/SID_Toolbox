clc; clear; close all hidden

%% Import Functions
addpath('../auxiliary_functions')

% Load & Save Results of the Surrogates for Significance and Non Linear
% Dynamics

% Significance
load('MIR_Significance_2CH_1.mat');
%%% KNN Estimation
MIR_Surr_KNN_qntl=quantile(KNN_MIR_Surr,0.95,3);
MIR_SURR_KNN_Results=KNN_MIR>MIR_Surr_KNN_qntl;
Significant_Surr_Results=sum(MIR_SURR_KNN_Results);

% Distribution MIR Values
x_plot = eps;
MIR_Values=KNN_MIR;

% Non Linear Dynamics
load('MIR_IAAFT_2CH_1.mat');
% KNN Estimation
MIR_Surr_qntl=quantile(MIR_Surr,0.95,3);

%% Compare Surr with Estimation
% KNN Estimation
MIR_SURR_Results=MIR_Est>MIR_Surr_qntl;
NonLinear_Surr_Results=sum(MIR_SURR_Results);

% Clear all variables execpt the results
clearvars -except Significant_Surr_Results MIR_Values NonLinear_Surr_Results x_plot

% Plot
fig=figure;
fig.WindowState='maximized';
subplot_tight(1,2,1,[0.16 0.08]);
errorbar(x_plot, median(MIR_Values), median(MIR_Values)-quantile(MIR_Values,0.05),quantile(MIR_Values,0.95)-median(MIR_Values),'LineWidth',2,'Color','k');
ylabel('I_{X;Y} [nats]','FontWeight','bold');
% ylim([-0.1 2.6]);
xlabel('$C$','FontWeight','bold','Interpreter','latex');
xticks(x_plot);
xlim([-0.01 0.51]);
ax=gca;
ax.FontSize=21;

subplot_tight(1,2,2,[0.16 0.08]);
Data=[Significant_Surr_Results' NonLinear_Surr_Results'];
bar(x_plot,Data);
ylim([0 115]);
xlabel('$C$','FontSize',25,'FontWeight','bold','Interpreter','latex');
xticks(x_plot);
xlim([-0.01 0.51]);
ylabel('%','FontSize',25,'FontWeight','bold');
yline(5,'-','LineWidth',2,'Color','r');
yticks(0:20:100);
% title('$$linear\ Model: Information\ Storage\ significance - removed\ nonlinearity$$','Interpreter','latex','FontSize',25)
legend({'Presence of Coupling' 'Presence of Nonlinearities'},'Box','off','Location','north','NumColumns',2);
ax=gca;
ax.FontSize=21;
AddLetters2Plots(fig,{'(a)','(b)'} ,'HShift',[0,0], 'VShift', 0, 'Direction', 'LeftRight','FontSize',20);

% Save Figure
exportgraphics(fig,'MIR_CH_KNN_errorbar_bar_new.eps','BackgroundColor','none');