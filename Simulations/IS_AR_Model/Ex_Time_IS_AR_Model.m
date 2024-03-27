%% Execution Time Study
clc; clear; close all hidden

%% Import Functions
addpath('../SID/functions')
addpath('../auxiliary_functions')
rng('default');

%% Parameters
signal_length=[100 500 1000 1500 2000]; % length of simulated time series

% signal_length=500;
f=0.3; % frequency of stochastic oscillation

M=1; %n. of time series
p=[2 4 6 8]; % maximum lag
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=1; %variance of innovation processes

num_signals = 100;
% num_surrogates = 100;

% embedding vector
m=p*ones(1,M);
tau=ones(1,M);

k=10; % n. of neighbors
rng('default');
r_arr = 0.5; % strength of stochastic oscillation

% Cells to store times
t_est=zeros(length(signal_length),length(m),num_signals);
t_sig=zeros(length(signal_length),length(m),num_signals);
t_iafft=zeros(length(signal_length),length(m),num_signals);

% Estimation Var Model
r = r_arr; % strength of stochastic oscillation
par.poles=([r f]); % Oscillation

%%% Theoretical VAR process
[Am,Su]=var_simulations(M,par); % parameters

for n=1:length(signal_length)
    for emb=1:length(m)
        for n_sim=1:num_signals
            VL=surr_SetLag(m(emb),tau);

            %% Calculate Information Storage - linear/nonparametric estimator
            Un = mvnrnd(zeros(1,M),Su,signal_length(n));
            Y = var_filter(Am,Un);

            % Information Storage - knn
            tic
            surr_ISknn(Y,VL,1,k,0);
            t_est(n, emb,n_sim) = toc;

            % Surrogates for Significance
            tic;
            surr_ISknn(Y,VL,1,k,1);
            t_sig(n, emb, n_sim) = toc;

            tic
            surr_ISknn(Y,VL,1,k,2);
            t_iafft(n,emb,n_sim)=toc;
        end
    end
end

%% Percentiles 
prctile_est = prctile(t_est,[5 95],3);
prctile_sig = prctile(t_sig,[5 95],3);
prctile_iafft = prctile(t_iafft,[5 95],3);

median_est=median(t_est,3);
median_sig=median(t_sig,3);
median_iafft=median(t_iafft,3);

save Execution_Times_IS_AR_Model

%% Quantiles plots
Colors_Plot=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];

fig=figure('WindowState','maximized');
subplot_tight(1, 3,1,[0.12 0.08]);
for i=1:length(m)
    errorbar(signal_length+i*15, median_est(:,i),median_est(:,i)-prctile_est(:,i,1),prctile_est(:,i,2)-median_est(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['q=' num2str(m(i))]); hold on;
    scatter(signal_length+i*15,median_est(:,i),[],Colors_Plot(i,:),'filled');
end
xlabel('Signal Length','FontWeight','bold');
xlim([60 2100]);
xticks(signal_length);
ylabel('Execution Time [sec]','FontWeight','bold');
title('Estimation Only','FontWeight','bold');
ax=gca;
ax.FontSize=18;

subplot_tight(1,3,2,[0.12 0.08]);
for i=1:length(m)
    errorbar(signal_length+i*15, median_sig(:,i),median_sig(:,i)-prctile_sig(:,i,1),prctile_sig(:,i,2)-median_sig(:,i) ,'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['q=' num2str(m(i))]); hold on;
    scatter(signal_length+i*15,median_sig(:,i),[],Colors_Plot(i,:),'filled');
end
xlabel('Signal Length','FontWeight','bold');
xlim([60 2100]);
xticks(signal_length);
title('Significance Testing','FontWeight','bold');
ax=gca;
ax.FontSize=18;

subplot_tight(1,3,3,[0.12 0.08]);
for i=1:length(m)
    errorbar(signal_length+i*15, median_iafft(:,i),median_iafft(:,i)-prctile_iafft(:,i,1),prctile_iafft(:,i,2)-median_iafft(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['q=' num2str(m(i))]); hold on;
    scatter(signal_length+i*15,median_iafft(:,i),[],Colors_Plot(i,:),'filled','HandleVisibility','off');
end
xlabel('Signal Length','FontWeight','bold');
xlim([60 2100]);
xticks(signal_length);
title('Non Linearities Testing','FontWeight','bold');
legend('Box','off');
ax=gca;
ax.FontSize=18;