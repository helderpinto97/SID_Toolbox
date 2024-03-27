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

k=10; % n. of neighbors

% Fix the Seed
rng('default');

%% Pre Allocate For Speed
knnSE = nan(num_sim, 1);
knnSE_surr = nan(num_sim, length(N), num_surrogates);

hw1 = waitbar(0,'Simulation Loop...');
hw2= waitbar(0,'Delta...');

for d=1:length(delta)
    for isim=1:num_sim

        % Generate Map
        Y=Hennon_Map_AdNoise(N,N_Transient,delta(d));

        %% Calculate Information Storage - nonparametric estimator

        % Information Storage - knn
        out=surr_ISknn(Y,VL,1,k,0);
        knn_Sy=out.Sy;
        knnSE(isim,d) = knn_Sy;

        %% Build Surrogates - destroy nonlinear correlation, preserver linear; Calculate Information Storage - linear/knn estimator

        parfor i=1:num_surrogates
            out=surr_ISknn(Y,VL,1,k,2);
            knn_Sy_surr=out.Sy;
            knnSE_surr(isim,d,i) = knn_Sy_surr;
        end

        waitbar(isim/num_sim,hw1);

    end
    waitbar(d/length(delta),hw2);
end

hw1.delete();
hw2.delete();

% KNN Significance
tmp_prctile_dyn = prctile(knnSE_surr,95,3);
knnSE_significance = knnSE > tmp_prctile_dyn;

% Percentages
per_KNN=sum(knnSE_significance);

%% Save Results
save IS_HH_IAFFT

%% Plot Results
Data=per_KNN;
fig=figure;
fig.WindowState='maximized';
bar(categorical(delta),Data,'LineWidth',1.5);
xlabel('$\delta$','Interpreter','latex','FontWeight','bold');
ylabel('\%','FontWeight','bold','Interpreter','latex');
ylim([0 101]);
yline(5,'LineWidth',1.5,'Color','red','FontSize',20);
ax=gca;
ax.FontSize=20;

