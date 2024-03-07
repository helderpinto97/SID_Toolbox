clc
clear
close all

%% Add path of functions
addpath('functions')

%% Load data
% Sample file from <XXX> dataset [REF]
file_name = 'CRa-01.prn'; % controlled breathing at 0.2 Hz
data=load(file_name);

series_name{1}='RR';
series_name{2}='SAP';
series_name{3}='RESP';

%% Plot series
figure(1);
for i=1:3 %3 series
    subplot(3,1,i); plot(data(:,i));
    title(['subject 01, CRa, ' series_name{i}]);
    xlim([1 length(data)])
end

%% Parameters

% Number of surrogates
num_surrogates = 1000;

% Number of neighbors
k=10;

% Embedding Matrix
M=1;
emb_dim=2;
m_emb=emb_dim*ones(1,M);
tau=ones(1,M);
VL=surr_SetLag(m_emb,tau);

%% Estimation Information Measures

% Estimation of IS for the RR
Y = data(:,1);

% Normalize data for the KNN estimator
Y=zscore(Y);
out=surr_SEknn(Y,VL,1,k);
IS=out.Sy;

% surrogate generation
IS_SignSurr = nan(num_surrogates,1);
IS_NLSurr = nan(num_surrogates,1);
for i=1:num_surrogates
    % Significance Surr
    out=surr_SEknn(Y,VL,1,k,1);
    IS_SignSurr(i)=out.Sy;

    % Nonlin Surr
    out=surr_SEknn(Y,VL,1,k,2);
    IS_NLSurr(i) = out.Sy;
end

%% Estimate 95th percentiles
SignSurr_95 = prctile(IS_SignSurr,95);
NLSurr_95 = prctile(IS_NLSurr,95);

%% Determine measure significance
fprintf("Information Storage estimation: %.2f\n", IS)
fprintf("Information Storage - measure significance: %s\n", mat2str(IS > SignSurr_95))
fprintf("Information Storage - nonlinearity presence: %s\n", mat2str(IS > NLSurr_95))
