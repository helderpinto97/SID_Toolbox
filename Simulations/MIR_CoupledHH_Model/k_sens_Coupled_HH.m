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
k=[5 10 15 20]; % Neighbours
V = [[ones(m,1);2*ones(m,1)], [(1:m)';(1:m)']]; % Embedding

%% 
M=2;
p=2;
pmax=20;
q=20;
tau=ones(1,M);

%%% Storage of Results
KNN_MIR=nan(length(eps),length(k),n_sim);
KNN_MIR_Surr=nan(n_sim,length(eps),num_surr);


% Waitbars
hw1 = waitbar(0,'C loop...');
hw2 = waitbar(0,'signal loop...');


for ik=1:length(k)
    for ieps=1:length(eps)
        for isim=1:n_sim
            [x,y] = Coupled_Henon(N, eps(ieps),rho);
            Y=[x y];
            data_Hennon_Hennon{isim,ieps}=Y;

            % KNN Estimation
            % Normalize Data
            Y=normalize(Y);

            % KNN Estimates
            out=surr_MIRknn(Y,V,1,2,k(ik),0);
            KNN_MIR(ieps,ik,isim)=out.I_XY2;

            % Surrogates For Significance
            % parfor nsurr=1:num_surr
            %     out=surr_MIRknn(Y,V,1,2,k,1);
            %     KNN_MIR_Surr(isim,ieps,nsurr)=out.I_XY2;
            % end
            waitbar(isim/n_sim,hw2);
        end
        waitbar(ieps/length(eps),hw1)
    end
end

hw1.delete();
hw2.delete();

%% Save the Results
save k_sens_MIR_Significance_2CH_1.mat

%% Percentiles 
prctile_est = prctile(KNN_MIR,[5 95], 3);

median_est=median(KNN_MIR,3);

save k_Sens_MIR_CHH

%% Quantiles plots
Colors_Plot=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];

fig=figure('WindowState','maximized');
% subplot_tight(1, 3,1,[0.12 0.08]);
for i=1:length(k)
    errorbar(eps, median_est(:,i),median_est(:,i)-prctile_est(:,i,1),prctile_est(:,i,2)-median_est(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['k=' num2str(k(i))]); hold on;
    scatter(eps,median_est(:,i),[],Colors_Plot(i,:),'filled','HandleVisibility','off');
end
% plot(eps, linSE_CE,'LineWidth',1.5,'Color','k', 'DisplayName','Theoretical Value');
xlabel('\epsilon','FontWeight','bold');
xticks(eps);
ylabel('MIR [nats]','FontWeight','bold');
title('K Sensivity -- Coupled Henon Henon Model','FontWeight','bold');
legend('Box','off');
ax=gca;
ax.FontSize=18;
