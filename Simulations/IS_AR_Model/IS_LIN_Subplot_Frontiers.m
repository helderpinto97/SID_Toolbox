%% Create the 2x1 subplot for Frontiers Papers
clc; clear all; close all hidden

%% Import Functions
addpath('../auxiliary_functions')

% Load & Save Results of the Surrogates for Significance and Non Linear
% Dynamics

% Significance
load('IS_Lin_Model_Significance');
Significant_Surr_Results=sum(knnSE_significance,2);

% Distribution of the IS Values - KNN
x = r_arr;
IS_Values=knnSE';

% Distribution Theoretical values
IS_Values_Th=linSE_CE;

% Non Linear Dynamics
load('IS_Lin_Model_IAAFT');
NonLinear_Surr_Results=sum(knnSE_significance,2);

% Clear all variables execpt the results
clearvars -except Significant_Surr_Results IS_Values NonLinear_Surr_Results x IS_Values_Th

% Plot
% HTML Gray
htmlGray = [128 128 128]/255;
% Figure
fig=figure;
fig.WindowState='maximized';
subplot_tight(1,2,1,[0.16 0.08]);
errorbar(x, median(IS_Values), median(IS_Values)-quantile(IS_Values,0.05,1),quantile(IS_Values,0.95,1)-median(IS_Values),'LineWidth',2,'Color','k','DisplayName','KNN');
hold on;
plot(x,IS_Values_Th,'LineWidth',2,'Color',htmlGray,'DisplayName','Theoretical Value');
hold off;
legend('Box','off','Location','best');
ylabel('S_X [nats]','FontWeight','bold');
xlabel('\rho_x','FontWeight','bold');
xticks(0:0.05:0.95);
xlim([-0.02 0.97]);
ax=gca;
ax.FontSize=21;

subplot_tight(1,2,2,[0.16 0.08]);
Data=[Significant_Surr_Results NonLinear_Surr_Results];
bar(x,Data);
ylim([0 110]);
xlabel('\rho_x','FontSize',25,'FontWeight','bold');
xticks(0:0.05:0.95);
xlim([-0.02 0.97]);
ylabel('%','FontSize',25,'FontWeight','bold');
yline(5,'-','LineWidth',2,'Color','r');
yticks(0:20:100);
% title('$$linear\ Model: Information\ Storage\ significance - removed\ nonlinearity$$','Interpreter','latex','FontSize',25)
legend({'Presence of Self-Dependencies' 'Presence of Nonlinearities'},'Box','off','Location','north','NumColumns',2);
ax=gca;
ax.FontSize=21;
AddLetters2Plots(fig,{'(a)','(b)'} ,'HShift',[0,0], 'VShift', 0, 'Direction', 'LeftRight','FontSize',18);

% Save Figure
exportgraphics(fig,'IS_LinearModel_KNN_errorbar_bar_new.eps','BackgroundColor','none');