%% Coupled Hennon-Hennon MAP - Non Linear Case
clc; clear; close all hidden

%% Import Functions
addpath('../SID/functions')
addpath('../auxiliary_functions')

%%% Number of Simulations
n_sim=100;

%%% Number of Surrogates
num_surr=100;

%%% Length of Series
N=1000;

%%% Coupling Parameter
eps = 0:0.05:0.5; % High coupling Parameter just before synchronization
rho=[0.3 0.3];

% Simulate Series
data_Hennon_Hennon=cell(length(eps),n_sim);

% Simulated Series Allocation
data_Hennon_Hennon=cell(length(eps),n_sim);

%%% KNN Estimation of MIR
%%% Parameters
m=2; % Embebbeding Dimension
k=10; % Neighbours
V = [[ones(m,1);2*ones(m,1)], [(1:m)';(1:m)']]; % Embedding

%% 
M=2;
p=2;
pmax=20;
q=20;
tau=ones(1,M);

%%% Storage of Results
KNN_MIR=nan(n_sim,length(eps));
KNN_MIR_Surr=nan(n_sim,length(eps),num_surr);


% Waitbars
hw1 = waitbar(0,'C loop...');
hw2 = waitbar(0,'signal loop...');

for ieps=1:length(eps)
    for isim=1:n_sim
        [x,y] = Coupled_Henon(N, eps(ieps),rho);
        Y=[x y];
        data_Hennon_Hennon{isim,ieps}=Y;

        % KNN Estimation
        % Normalize Data
        Y=normalize(Y);

        % KNN Estimates 
        out=surr_MIRknn(Y,V,1,2,k,0);
        KNN_MIR(isim,ieps)=out.I_XY2;

        % Surrogates For Significance 
        parfor nsurr=1:num_surr
            out=surr_MIRknn(Y,V,1,2,k,1);
            KNN_MIR_Surr(isim,ieps,nsurr)=out.I_XY2;
        end
          waitbar(isim/n_sim,hw2);
    end
    waitbar(ieps/length(eps),hw1)
end

hw1.delete();
hw2.delete();

%% Save the Results
save MIR_Significance_2CH_1.mat

%% NON Linear Model Plots

%%% KNN Estimation
MIR_Surr_KNN_qntl=quantile(KNN_MIR_Surr,0.95,3);
MIR_SURR_KNN_Results=KNN_MIR>MIR_Surr_KNN_qntl;
MIR_SURR_KNN_Results_Sum=sum(MIR_SURR_KNN_Results);

% Plot Results 
fig=figure;
fig.WindowState='maximized';
bar(eps,MIR_SURR_KNN_Results_Sum);
xlim([-0.05 0.55]);
xlabel('$C$','FontSize',20,'FontWeight','bold','Interpreter','latex');
xticks(eps);
ylabel('%','FontSize',20,'FontWeight','bold');
% title({'KNN Estimation', 'Percentage of Significant Measures'});
yline(5,'-','LineWidth',2,'Color','r');
ax = gca;
ax.FontSize=20;

% Mean Profiles
mean_KNN_MIR=mean(KNN_MIR);

std_KNN_MIR=3*std(KNN_MIR);

fig2=figure;
fig2.WindowState='maximized';
errorbar(eps,mean_KNN_MIR,std_KNN_MIR,'LineWidth',1.5,'Color','r','DisplayName','KNN Estimator'); hold on;
area(eps, mean(MIR_Surr_KNN_qntl,1),'FaceColor','r','EdgeColor','none','FaceAlpha',0.4,'DisplayName','95% Percentiles Surrogates KNN');
legend('Box','off','Location','best');
xlabel('$\epsilon$','Interpreter','latex');
ylabel('$I_{X;Y}$','Interpreter','latex');
ax=gca;
ax.FontSize=20;

%% Save Figures
% exportgraphics(fig1,'MIR_NLinearModel_Significance.png','BackgroundColor','none');
% exportgraphics(fig2,'MIR_NLinearModel_MeanProf.png','BackgroundColor','none');

 