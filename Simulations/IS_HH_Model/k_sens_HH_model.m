clc; clear; close all hidden

%% Import Functions
addpath('../SID/functions')
addpath('../auxiliary_functions')

%% Parameters
delta=0:0.2:3;
N=1000;
N_Transient=1000;

num_sim = 100;
num_surrogates = 100;

% embedding vector
M=1; %n. of time series
p=2; % maximum lag
pmax=20;
m=p*ones(1,M);
tau=ones(1,M);
VL=surr_SetLag(m,tau);

k=[5 10 15 20]; % n. of neighbors

% Fix the Seed
rng('default');

%% Pre Allocate For Speed
knnSE = nan(length(delta),length(k),num_sim);
% knnSE_surr = nan(num_sim, length(N), num_surrogates);

hw1 = waitbar(0,'Simulation Loop...');
hw2= waitbar(0,'Delta...');
hw3= waitbar(0,'K...');

for ik=1:length(k)
    for d=1:length(delta)
        for isim=1:num_sim

            % Generate Map
            Y=Hennon_Map_AdNoise(N,N_Transient,delta(d));

            %% Calculate Information Storage - linear/nonparametric estimator

            % Information Storage - knn
            out=surr_ISknn(Y,VL,1,k(ik),0);
            knn_Sy=out.Sy;
            knnSE(d,ik,isim) = knn_Sy;

            %% Build Surrogates - destroy nonlinear correlation, preserver linear; Calculate Information Storage - linear/knn estimator

            % parfor i=1:num_surrogates
            %     out=surr_ISknn(Y,VL,1,k,1);
            %     knn_Sy_surr=out.Sy;
            %     knnSE_surr(isim,d,i) = knn_Sy_surr;
            % end

            waitbar(isim/num_sim,hw1);

        end
        waitbar(d/length(delta),hw2);
    end
end

hw1.delete();
hw2.delete();
hw3.delete();

% Percentiles
prctile_est = prctile(knnSE,[5 95],3);
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
    errorbar(delta, median_est(:,i),median_est(:,i)-prctile_est(:,i,1),prctile_est(:,i,2)-median_est(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['k=' num2str(k(i))]); hold on;
    scatter(delta,median_est(:,i),[],Colors_Plot(i,:),'filled','HandleVisibility','off');
end
% plot(r_arr, linSE_CE,'LineWidth',1.5,'Color','k', 'DisplayName','Theoretical Value');
xlabel('\delta','FontWeight','bold');
xticks(delta);
ylabel('IS [nats]','FontWeight','bold');
title('K Sensivity -- Modify Henon Henon Model','FontWeight','bold');
legend('Box','off');
ax=gca;
ax.FontSize=18;

