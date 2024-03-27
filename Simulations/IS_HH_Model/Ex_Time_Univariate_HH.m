%% Execution Times for the Univariate
clc; clear; close all hidden

%% Import Functions
addpath('../SID/functions')
addpath('../auxiliary_functions')

%% Parameters
% delta=0:0.2:3;

% For this study define only a value of delta
delta=1;

N=[100 500 1000 1500 2000];
N_Transient=1000;

num_sim = 100;
num_surrogates = 100;

% embedding vector
M=1; %n. of time series
p=[2 4 6 8]; % maximum lag
pmax=20;
m=p*ones(1,M);
tau=ones(1,M);

k=10; % n. of neighbors

% Fix the Seed
rng('default');


%% Pre Allocate For Speed
t_est=nan(length(N),length(m),num_sim);
t_sig=nan(length(N),length(m),num_sim);
t_iafft=nan(length(N),length(m),num_sim);

hw1 = waitbar(0,'Simulation Loop...');
hw2= waitbar(0,'Embedding Dimension...');
hw3= waitbar(0,'Length Signal...');

for n=1:length(N)
    for emb=1:length(m)

        % Define Embedding matrix
        VL=surr_SetLag(m(emb),tau);
        
        for isim=1:num_sim

            % Generate Map
            Y=Hennon_Map_AdNoise(N(n),N_Transient,delta);

            %% Estimation only
            tic;
            out=surr_ISknn(Y,VL,1,k,0);
            t_est(n,emb,isim)=toc;

            %% Surrogates for Significance
            tic;
            out=surr_ISknn(Y,VL,1,k,1);
            t_sig(n,emb,isim)=toc;

            %% Surrogates for NonLinearities
            tic;
            out=surr_ISknn(Y,VL,1,k,2);
            t_iafft(n,emb,isim)=toc;


            waitbar(isim/num_sim,hw1);

        end
        waitbar(emb/length(m),hw2);
    end
    waitbar(n/length(N),hw3);
end

hw1.delete();
hw2.delete();
hw3.delete();

%% Plots Results

%% Percentiles 
prctile_est = prctile(t_est,[5 95],3);
prctile_sig = prctile(t_sig,[5 95],3);
prctile_iafft = prctile(t_iafft,[5 95],3);

median_est=median(t_est,3);
median_sig=median(t_sig,3);
median_iafft=median(t_iafft,3);

save Execution_Times_IS_HH_Model

%% Quantiles plots
Colors_Plot=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];

fig=figure('WindowState','maximized');
subplot_tight(1, 3,1,[0.12 0.08]);
for i=1:length(m)
    errorbar(N+i*15, median_est(:,i),median_est(:,i)-prctile_est(:,i,1),prctile_est(:,i,2)-median_est(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['q=' num2str(m(i))]); hold on;
    scatter(N+i*15,median_est(:,i),[],Colors_Plot(i,:),'filled');
end
xlabel('Signal Length','FontWeight','bold');
xlim([60 2100]);
xticks(N);
ylabel('Execution Time [sec]','FontWeight','bold');
title('Estimation Only','FontWeight','bold');
ax=gca;
ax.FontSize=18;

subplot_tight(1,3,2,[0.12 0.08]);
for i=1:length(m)
    errorbar(N+i*15, median_sig(:,i),median_sig(:,i)-prctile_sig(:,i,1),prctile_sig(:,i,2)-median_sig(:,i) ,'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['q=' num2str(m(i))]); hold on;
    scatter(N+i*15,median_sig(:,i),[],Colors_Plot(i,:),'filled');
end
xlabel('Signal Length','FontWeight','bold');
xlim([60 2100]);
xticks(N);
title('Significance Testing','FontWeight','bold');
ax=gca;
ax.FontSize=18;

subplot_tight(1,3,3,[0.12 0.08]);
for i=1:length(m)
    errorbar(N+i*15, median_iafft(:,i),median_iafft(:,i)-prctile_iafft(:,i,1),prctile_iafft(:,i,2)-median_iafft(:,i),'Color',Colors_Plot(i,:),'LineWidth',2,'DisplayName',['q=' num2str(m(i))]); hold on;
    scatter(N+i*15,median_iafft(:,i),[],Colors_Plot(i,:),'filled','HandleVisibility','off');
end
xlabel('Signal Length','FontWeight','bold');
xlim([60 2100]);
xticks(N);
title('Non Linearities Testing','FontWeight','bold');
legend('Box','off');
sgtitle('Univariate Hennon Map ','FontWeight','bold','FontSize',18);
ax=gca;
ax.FontSize=18;