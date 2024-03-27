%% Hypothesis 02a - Underlying dynamic is linear in 2AR model (bivariate case)
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
m=p.*ones(1,M);
tau=ones(1,M);
VL=surr_SetLag(m,tau);

k=10; % n. of neighbors

rng('default');

%% Experiment 01 - variable coupling strength of one processes
disp('Experiment 01 - variable coupling strength of one processes')
tic

N=1000; % length of simulated time series
C_arr = 0:0.05:1; % variable coupling strength of one processes

knnMIR_varC = nan(length(C_arr), num_signals, 1);
knnMIR_surr_varC = nan(length(C_arr), num_signals, num_surrogates);

linMIR_theoretical_varC = nan(length(C_arr), 1);

% Build Simulation
hw1 = waitbar(0,'C loop...');
hw2 = waitbar(0,'signal loop...');

COH_org = {};
COH_surr = {};
q=10;
for C_idx = 1:length(C_arr) 
    C = C_arr(C_idx); % strength of stochastic oscillation

    par.coup=[1 2 2 C]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
    %%% Theoretical VAR process
    [Am,Su]=var_simulations(M,par); % parameters

    ret = surr_mir_th(Am,Su,q,2,1);
    linMIR_theoretical_varC(C_idx) = ret.I_XY2;
    
    for sig_idx = 1:num_signals
        
        % Estimation on a realization of the simulation
        Un = mvnrnd(zeros(1,M),Su,N);
        Y = var_filter(Am,Un); % realization
        Y = zscore(Y);
        
        % Calculate Mutual Information Rate - linear/nonparametric estimator
        
        % Mutual Information Rate - knn  
        out=surr_MIRknn(Y,VL,1,2,k,0);
        knnMIR_varC(C_idx, sig_idx) = out.I_XY2;

        % Build Surrogates - remove nonlinearities; Calculate MIR - linear/knn estimator
        parfor i = 1:num_surrogates
           % Surrogates Generation - IAAFT
            outs=surr_MIRknn(Y,VL,1,2,k,2);

            knnMIR_surr_varC(C_idx, sig_idx, i) = outs.I_XY2;
        end
        clear out
        waitbar(sig_idx/num_signals,hw2);
    end

    waitbar(C_idx/length(C_arr),hw1);
end

hw1.delete;
hw2.delete;

disp('done')
toc
save MIR_VAR_IAFFT

%% plots
tmp_prctile_knn = prctile(knnMIR_surr_varC,95,3);
knnSE_significance = knnMIR_varC > tmp_prctile_knn;

x = C_arr;
figure('WindowState', 'maximized')
area(x, mean(tmp_prctile_knn,2),'FaceColor','r','EdgeColor','none','FaceAlpha',0.4,'DisplayName','95% Percentile Surrogates KNN Estimator');
hold on
errorbar(x, mean(knnMIR_varC,2)', 2*std(knnMIR_varC,[],2)','LineWidth',2,'Color','r','DisplayName','KNN Estimator')
plot(x, linMIR_theoretical_varC,'LineWidth',2,'Color','k','Marker','x','MarkerSize',10,'DisplayName','Theoretical value')
xlabel('$$C$$','Interpreter','latex','FontSize',25)
ylabel('$$MIR\ [nats]$$','Interpreter','latex','FontSize',25)
%title('$$2AR\ Model: Mutual\ Information\ Rate\ - removed\ nonlinearities$$','Interpreter','latex','FontSize',25)
legend('Location','bestoutside');
ax=gca;
ax.FontSize=20;
% exportgraphics(gcf,[out_fig_path, 'MIR_2AR_nonliniarity_C_errorbar_V2.png'],'BackgroundColor','none');

% Only KNN Estimator
figure('WindowState', 'maximized')
bplot = bar(x,mean(knnSE_significance,2)*100);
% bplot(1).FaceColor = 'b';
% bplot(2).FaceColor = 'r';
% set(bplot,{'DisplayName'},{'Linear Estimator'; 'KNN Estimator'})
ylim([0 100])
xlabel('$$C$$','Interpreter','latex','FontSize',25,'FontWeight','bold');
ylabel('$$\%$$','Interpreter','latex','FontSize',25,'FontWeight','bold');
xticks(x);
%title('$$2AR\ Model: Mutual\ Information\ Rate\ significance - removed\ nonlinearities$$','Interpreter','latex','FontSize',25)
yline(5,'Linewidth',2,'Color','r', 'DisplayName','5%','HandleVisibility','off');
% legend('Location','bestoutside');
ax=gca;
ax.FontSize=20;