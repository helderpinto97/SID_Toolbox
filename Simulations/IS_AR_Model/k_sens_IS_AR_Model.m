%% Sensivity to number of neighbors k
clc; clear; close all hidden

%% Import Functions
addpath('../SID/functions')
addpath('../auxiliary_functions')
rng('default');

%% Parameters
signal_length=1000; % length of simulated time series

% signal_length=500;
f=0.3; % frequency of stochastic oscillation

M=1; %n. of time series
p=2; % maximum lag
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=1; %variance of innovation processes

num_signals = 100;
% num_surrogates = 100;

% embedding vector
m=p*ones(1,M);
tau=ones(1,M);
VL=surr_SetLag(m,tau);

k=[5 10 15 20]; % n. of neighbors
rng('default');
r_arr = 0:0.05:0.95; % strength of stochastic oscillation
num_surrogates=100;

%% variables
linSE_CE = nan(length(r_arr), 1);

knnSE = nan(length(r_arr), length(k), num_signals);
knnSE_surr = nan(length(r_arr), length(k),num_signals, num_surrogates);

%% Build Simulation

hw1 = waitbar(0,'r loop...');
hw2 = waitbar(0,'signal loop...');
hw3= waitbar(0,'Number of neighbors');


for ik=1:length(k)
    for r_idx = 1:length(r_arr)
        r = r_arr(r_idx); % strength of stochastic oscillation
        par.poles=([r f]); % Oscillation

        %%% Theoretical VAR process
        [Am,Su]=var_simulations(M,par); % parameters

        ret = surr_CElinVAR1(Am,Su,2);
        linSE_CE(r_idx) = ret.Hy-ret.Hy_y;

        for sig_idx = 1:num_signals
            %% Calculate Information Storage - linear/nonparametric estimator
            Un = mvnrnd(zeros(1,M),Su,signal_length);
            Y = var_filter(Am,Un);

            % Information Storage - knn
            out=surr_ISknn(Y,VL,1,k(ik),0);
            knnSE(r_idx, ik,sig_idx) = out.Sy;

            %% Build Surrogates - destroy the dynamics; Calculate Information Storage - knn estimator
            % parfor i = 1:num_surrogates
            %     out=surr_ISknn(Y,VL,1,k(ik),1);
            %     knnSE_surr(r_idx, ik,sig_idx,i) = out.Sy;
            % end

            waitbar(sig_idx/num_signals,hw2);
        end

        waitbar(r_idx/length(r_arr),hw1);
    end
    waitbar(ik/length(k),hw3);
end

hw1.delete;
hw2.delete;
hw3.delete;


%% Percentiles 
prctile_est = prctile(knnSE,[5 95], 3);

median_est=median(knnSE,3);

save k_Sens_IS_AR_Model

%% Quantiles plots
Colors_Plot=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];

fig=figure('WindowState','maximized');
% subplot_tight(1, 3,1,[0.12 0.08]);
for i=1:length(k)
    errorbar(r_arr, median_est(:,i),median_est(:,i)-prctile_est(:,i,1),prctile_est(:,i,2)-median_est(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['k=' num2str(k(i))]); hold on;
    scatter(r_arr,median_est(:,i),[],Colors_Plot(i,:),'filled','HandleVisibility','off');
end
plot(r_arr, linSE_CE,'LineWidth',1.5,'Color','k', 'DisplayName','Theoretical Value');
xlabel('\rho_x','FontWeight','bold');
xticks(r_arr);
ylabel('IS [nats]','FontWeight','bold');
title('K Sensivity -- AR Model','FontWeight','bold');
legend('Box','off');
ax=gca;
ax.FontSize=18;

