clc
clear
close all

%% Path to Functions
addpath('functions')

%% Load data

% Files from Faes, L., Nollo, G., and Porta, A. (2011). Information domain approach to the investigation of cardio-586
% vascular, cardio-pulmonary, and vasculo-pulmonary causal couplings. Frontiers in Physiology 2 NOV.587 doi:10.3389/fphys.2011.00080

file_name = 'CRa-01.prn'; %controlled breathing at 0.2 Hz
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
M=2;
emb_dim=2;
m_emb=emb_dim*ones(1,M);
tau=ones(1,M);
VL=surr_SetLag(m_emb,tau);

%% Estimation Information Measures

% Estimation of MIR between RR & RESP
Y = data(:,[1,3]);

% Normalize data for the KNN estimator
Y=zscore(Y);
out=surr_MIRknn(Y,VL,1,2,k);
MIR=out.I_XY2;

% surrogate generation
MIR_SignSurr = nan(num_surrogates,1);
MIR_NLSurr = nan(num_surrogates,1);
for is=1:num_surrogates
    % Significance Surr
    out=surr_MIRknn(Y,VL,1,2,k,1);
    MIR_SignSurr(is)=out.I_XY2;

    % Nonlin Surr
    out=surr_MIRknn(Y,VL,1,2,k,2);
    MIR_NLSurr(is) = out.I_XY2;
end

%% Estimate 95th percentiles
SignSurr_95 = prctile(MIR_SignSurr,95);
NLSurr_95 = prctile(MIR_NLSurr,95);

%% Determine measure significance

fprintf("Mutual Information Rate estimation: %.2f\n", MIR)
fprintf("Mutual Information Rate - measure significance: %s\n", mat2str(MIR > SignSurr_95))
fprintf("Mutual Information Rate - nonlinearity presence: %s\n", mat2str(MIR > NLSurr_95))
