clc; clear; close all hidden

%% Import Functions
addpath('../SID/functions')
addpath('../auxiliary_functions')

%% parameters
M=2;
p=3;
pmax = 20;
par.poles=([0.3 0.3; 0.3 0.1]); % Oscillations
par.Su=[1 1]; %variance of innovation processes

nit = 50; % iaafft iteration number

num_signals = 100;
num_surrogates = 100;

fs=1; % sampling frequency
nfft=1001; % number of points on frequency axis (total)

% embedding vector
m=2;
tau=ones(1,M);
VL=surr_SetLag(m,tau);
k=[5 10 15 20]; % n. of neighbors


rng('default');

%% Experiment 01 - variable coupling strength of one processes
disp('Experiment 01 - variable coupling strength of one processes')
tic

N=1000; % length of simulated time series
C_arr = 0:0.05:1; % variable coupling strength of one processes

% saved variables
knnMIR_varC = nan(length(C_arr), length(k),num_signals);

linMIR_theoretical_varC = nan(length(C_arr), 1);

% Build Simulation
hw1 = waitbar(0,'C loop...');
hw2 = waitbar(0,'signal loop...');
hw3= waitbar(0,'K - Neighbors Number...');
q=20;

for ik=1:length(k)
    for C_idx = 1:length(C_arr)
        C = C_arr(C_idx); % strength of stochastic oscillation

        par.coup=[1 2 2 C]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c
        %%% Theoretical VAR process
        [Am,Su]=var_simulations(M,par); % parameters

        ret =surr_mir_th(Am,Su,q,2,1);
        linMIR_theoretical_varC(C_idx) = ret.I_XY2;

        for sig_idx = 1:num_signals

            % Estimation on a realization of the simulation
            Un = mvnrnd(zeros(1,M),Su,N);
            Y = var_filter(Am,Un); % realization
            Y = zscore(Y);

            % Calculate Mutual Information Rate - nonparametric estimator

            % Mutual Information Rate - knn
            out=surr_MIRknn(Y,VL,1,2,k(ik),0);
            knnMIR_varC(C_idx,ik,sig_idx) = out.I_XY2;

            % Build Surrogates - remove nonlinearities; Calculate MIR - linear/knn estimator
            % parfor i = 1:num_surrogates
            %     % Surrogates Generation - Significance
            %     % KNN Estimation
            %     out=surr_MIRknn(Y,VL,1,2,k,1);
            %     knnMIR_surr_varC(C_idx,sig_idx,i) = out.I_XY2;
            % end
            waitbar(sig_idx/num_signals,hw2);
        end
        waitbar(C_idx/length(C_arr),hw1);
    end
    waitbar(ik/length(k),hw3);
end

hw1.delete;
hw2.delete;
hw3.delete;

disp('done')
toc
save k_sens_MIR_VAR_Model

prctile_est = prctile(knnMIR_varC,[5 95], 3);
median_est=median(knnMIR_varC,3);

%% Quantiles plots
Colors_Plot=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];

fig=figure('WindowState','maximized');
% subplot_tight(1, 3,1,[0.12 0.08]);
for i=1:length(k)
    errorbar(C_arr, median_est(:,i),median_est(:,i)-prctile_est(:,i,1),prctile_est(:,i,2)-median_est(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['k=' num2str(k(i))]); hold on;
    scatter(C_arr,median_est(:,i),[],Colors_Plot(i,:),'filled','HandleVisibility','off');
end
plot(C_arr, linMIR_theoretical_varC,'LineWidth',1.5,'Color','k', 'DisplayName','Theoretical Value');
xlabel('\rho_x','FontWeight','bold');
xticks(C_arr);
ylabel('MIR [nats]','FontWeight','bold');
title('K Sensivity -- VAR Model','FontWeight','bold');
legend('Box','off');
ax=gca;
ax.FontSize=18;