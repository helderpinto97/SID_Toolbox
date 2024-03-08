clc; clear all; close all

%% Import Functions
addpath('..\functions')


%% parameters
signal_length=1000; % length of simulated time series
f=0.3; % frequency of stochastic oscillation

M=1; %n. of time series
p=2; % maximum lag
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=1; %variance of innovation processes

num_signals = 100;
num_surrogates = 100;

% embedding vector
m=p*ones(1,M);
tau=ones(1,M);
VL=surr_SetLag(m,tau);

k=10; % n. of neighbors

rng('default');

r_arr = 0:0.05:0.95; % strength of stochastic oscillation

%% variables

linSE = nan(length(r_arr), num_signals, 1);
linSE_surr = nan(length(r_arr), num_signals, num_surrogates);

knnSE = nan(length(r_arr), num_signals, 1);
knnSE_surr = nan(length(r_arr), num_signals, num_surrogates);

%% Build Simulation

hw1 = waitbar(0,'r loop...');
hw2 = waitbar(0,'signal loop...');

for r_idx = 1:length(r_arr) 
    r = r_arr(r_idx); % strength of stochastic oscillation

    par.poles=([r f]); % Oscillation
    %%% Theoretical VAR process
    [Am,Su]=var_simulations(M,par); % parameters
    
    for sig_idx = 1:num_signals
        
        % Estimation on a realization of the simulation
        Un = mvnrnd(zeros(1,M),Su,signal_length);
        Y = var_filter(Am,Un); % realization
        
        %% Calculate Information Storage - nonparametric estimator
        
        % Information Storage - knn
        out=surr_ISknn(Y,VL,1,k,0);     
        knnSE(r_idx, sig_idx) = out.Sy;
        
        %% Build Surrogates - destroy nonlinear correlation, preserver linear; Calculate Information Storage - linear/knn estimator
        
        parfor i = 1:num_surrogates
            out=surr_ISknn(Y,VL,1,k,2);
            knnSE_surr(r_idx, sig_idx, i) = out.Sy;
        end
    
        waitbar(sig_idx/num_signals,hw2);
    end

    waitbar(r_idx/length(r_arr),hw1);
end

hw1.delete;
hw2.delete;
%%

tmp_prctile = prctile(linSE_surr,95,3);
linSE_significance = linSE > tmp_prctile;

tmp_prctile = prctile(knnSE_surr,95,3);
knnSE_significance = knnSE > tmp_prctile;

save IS_Lin_Model_IAAFT

%%
x = r_arr;
figure('WindowState', 'maximized')
area(x, mean(tmp_prctile,2),'FaceColor','g','EdgeColor','none')
hold on
errorbar(x, mean(knnSE,2), 3*std(knnSE,[],2),'LineWidth',2,'Color','r');
legend({'95% Percentile Surrogates' 'IS KNN'},'Box','off');
xlabel('$$\rho$$','Interpreter','latex','FontSize',25)
ylabel('$$IS\ [nats]$$','Interpreter','latex','FontSize',25)
% title('$$linear\ Model: Information\ Storage\ - removed\ nonlinearity$$','Interpreter','latex','FontSize',25)
ax=gca;
ax.FontSize=20;


figure('WindowState', 'maximized');
bar(x, sum(knnSE_significance,2));
ylim([0 100]);
xlabel('$$\rho$$','Interpreter','latex','FontSize',25);
xticks(x);
ylabel('%','Interpreter','tex','FontSize',25);
yline(5,'-','LineWidth',2,'Color','r');
ax=gca;
ax.FontSize=20;