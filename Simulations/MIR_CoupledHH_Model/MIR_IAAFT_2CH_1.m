%% Verify if the same results using the Hanning Window for the non linear model
clear all; close all; clc;


addpath('../../functions')

%%% Number of Surrogates
num_surr=100;

%%% Number of Simulations
n_sim=100;

%%% Coupling Parameter
eps = 0:0.05:0.5; % High coupling Parameter just before synchronization
rho=[0.3 0.3];

%%% Sample Size
N=1000;

%%% Max VAR order
p_max=20;

%%% Embedding Dimension
emb=2;

% Simulate Series
data_Hennon_Hennon=cell(n_sim,length(eps));

for ieps=1:length(eps)
    for isim=1:n_sim
        [x,y] = Coupled_Henon(N, eps(ieps),rho);
        Y=[x y];
        data_Hennon_Hennon{isim,ieps}=Y;
    end
end

%% KNN Estimation of MIR
%%% Parameters
M=2;
m=2; % Embebbeding Dimension
k=10; % Neighbours
V = [[ones(m,1);2*ones(m,1)], [(1:m)';(1:m)']]; % Embedding
tau=ones(1,M);
q=20;

MIR_Est=nan(n_sim,length(eps));

for ieps=1:length(eps)
    for isim=1:length(data_Hennon_Hennon)
        Yo=data_Hennon_Hennon{isim,ieps};

        % Normalize Data
        Yo=normalize(Yo);

        % Estimates
        out=surr_MIRknn(Yo,V,1,2,k,0);
        MIR_Est(isim,ieps)=out.I_XY2;
    end
end

%% Generate Surrogates
nit=50;
iX1=1;
iX2=2;

%%% Allocate Surrogates Data
MIR_Surr=nan(n_sim,length(eps),num_surr);
MIR_Linear_Surr=nan(n_sim,length(eps),num_surr);

% Build Surrogates
hw1 = waitbar(0,'C loop...');
hw2 = waitbar(0,'signal loop...');

for ieps=1:length(eps)
    for isim=1:n_sim
        DATA=data_Hennon_Hennon{isim,ieps};
        parfor n_surr=1:num_surr
            % Normalize for the KNN Estimation
            Y=normalize(DATA);
            % MIR - KNN Estimation
            out=surr_MIRknn(Y,V,1,2,k,2);
            MIR_Surr(isim,ieps,n_surr)=out.I_XY2;
        end
        waitbar(isim/n_sim,hw2);
    end
    waitbar(ieps/length(eps),hw1);
end

hw1.delete;
hw2.delete;

save MIR_IAAFT_2CH_1

%% Estimate Quantiles of Surrogates

% KNN Estimation
MIR_Surr_qntl=quantile(MIR_Surr,0.95,3);


%% Compare Surr with Estimation
% KNN Estimation
MIR_SURR_Results=MIR_Est>MIR_Surr_qntl;
MIR_SURR_Results_Sum=sum(MIR_SURR_Results);


% Plot Results KNN Estimator
fig=figure;
fig.WindowState='maximized';
bar(eps,MIR_SURR_Results_Sum);
xlim([-0.05 0.55]);
xlabel('$C$','FontSize',18,'FontWeight','bold','Interpreter','latex');
xticks(eps);
ylabel('%','FontSize',18,'FontWeight','bold');
yline(5,'LineWidth',2,'Color','red','HandleVisibility','off');
% title({'KNN Estimation', 'Percentage of Significant Measures'});
ax = gca;
ax.FontSize=20;

mean_KNN_MIR=mean(MIR_Est);
std_KNN_MIR=3*std(MIR_Est);

fig2=figure;
fig2.WindowState='maximized';
errorbar(eps,mean_KNN_MIR,std_KNN_MIR,'LineWidth',1.5,'Color','b','DisplayName','KNN Estimator'); hold on;
area(eps, mean(MIR_Surr_qntl,1),'FaceColor','b','EdgeColor','none','FaceAlpha',0.4,'DisplayName','95% Percentile Surrogates KNN');
legend('Box','off','Location','best');
xlabel('$\epsilon$','Interpreter','latex');
ylabel('$I_{X;Y}$','Interpreter','latex');
ax=gca;
ax.FontSize=20;