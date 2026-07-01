%%% A4_SUB_FDT_ANALYZE_DATA_1  Reduces FDT matrices to per-subject scalar metrics.
% Author: Juan Manuel Monti
%
% Loads partial_Sleep_COND_*.mat from a3 and integrates FDT observables over (t,s)
% via trapz to obtain one scalar per subject (averaged over parcels first).
%
% Reduction steps:
%   1. Average |XFDR-1|, |dVFDT|, |iVFDT| over parcels  -> [NSUBxTupxTup]
%   2. Integrate over t (trapz) for each fixed s          -> [NSUBxTup]
%   3. Integrate over s (trapz)                           -> [NSUBx1]
%
% Output metrics (one scalar per subject):
%   MX_sub   |XFDR - 1| integrated over (t,s)
%   Md_sub   |dVFDT|    integrated over (t,s)
%   Mi_sub   |iVFDT|    integrated over (t,s)
%   xMX_sub, xMd_sub, xMi_sub  same using alternative response (xR)
%
% Requires: partial_Sleep_COND_*.mat
% Saves:    metrics_SUB_Sleep_COND_<C>_NSUB_<N>_NSIM_<S>_DFILT_<D>_SIGMA_<SG>_LIN_<L>_HINIT_<H>_FREQSUB_<F>.mat

clearvars

Nsimulations = [100 250 750 1000 2500 7500 10000]
for NSUBSIM = Nsimulations
NSUBSIM

%%% Choose dataset and other parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SLEEP DATASET
% /Data/Sleep (15 subjects - 90 parcels)

%%% SET CONDITION
% FOR SLEEP DATASET:
% 1 --> Awake
% 2 --> Sleep
for CONDITION = [1 2];

%%% Use filtered or non-filtered data (ts, force, noise) to calculate C, A and R
%   DATAFILTER = 1 --> use filtered data
%   DATAFILTER = 0 --> use non-filtered data
DATAFILTER = 0;

%%% Set value of standard deviation of noise
if CONDITION == 1
    SETSIGMA = 0.12;
elseif CONDITION == 2
    SETSIGMA = 0.06;
end

%%% Linearization: Linear (LIN=1 NOT YET IMPLMENTED) or Cubic (LIN=0) terms in Hopf
LIN = 0;

%%% Hopf initialization
%   HOPFINT = 1 --> initializes for each individual simulation
%   HOPFINT = 0 --> initializes ONLY for the first simulation (of each subject)
%                   and all subsequent use the same z0
HOPFINIT = 0;

%%% Hopf Frequencies
%   FREQSUB = 1 --> different frequencies for each subject
%   FREQSUB = 0 --> frequencies for the mean Power Spectra
FREQSUB = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0. LOAD EMPIRIC DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structural Connectivity and Time series
load('./Data/Sleep/DataSleep_W_N3.mat')
% Nodes Frequencies and FC empiric, the output of Compute_Hopf_freq_AAL90
if CONDITION == 1
    load('./Data/Sleep/hopf_freq_AAL90_COND_1_W.mat') % the output of Compute_Hopf_freq_AAL90
    TS_X = TS_W;
elseif  CONDITION == 2
    load('./Data/Sleep/hopf_freq_AAL90_COND_2_N3.mat') % the output of Compute_Hopf_freq_AAL90
    TS_X = TS_N3;
end
%%% SC
C = SC/max(max(SC))*0.2;
%%% Diagonal of C equal zero
C = C.*~eye(size(C));

%%% Parcels
NPARCELS = length(C);
disp(['NPARCELS = ', num2str(NPARCELS)])

%%% Subjects in empiric data
NSUB = size(TS_X,2);
disp(['NSUB = ', num2str(NSUB)])
% Sets NSUBSIM = NSUB
if NSUBSIM == 0
    NSUBSIM = NSUB;
end
disp(['NSUBSIM = ', num2str(NSUBSIM)])

%%% Tmax & TR
Tmax = zeros(1,NSUB);
for i = 1:NSUB; Tmax(i) = size(TS_X{1,i},2); end
Tmax = min(min(Tmax));
% Tmax = 198;
disp(['Tmax = ', num2str(Tmax)])
Tinf = 10;
Tsup = Tmax - Tinf;
Tup = Tsup - Tinf + 1;
disp(['Tinf = ', num2str(Tinf)])
disp(['Tsup = ', num2str(Tsup)])
% Repetition time
TR = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. LOAD SIMULATED DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = sprintf('partial_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d.mat',...
                              CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
if isfile(file)
    load(file)
    disp('! Load simulated FDT data')
else
    disp('! No simulated data --> run FDT')
    return
end

%%% For each element I calculate the METRICS %%%
MX_all = abs(XFDRsub - 1);
Md_all = abs(dVFDTsub);
Mi_all = abs(iVFDTsub);
xMX_all = abs(xXFDRsub - 1);
xMd_all = abs(xdVFDTsub);
xMi_all = abs(xiVFDTsub);

%%% Average over parcels
% Author: Juan Manuel Monti
MX_sub_tt_ss = squeeze(mean(MX_all,2));
Md_sub_tt_ss = squeeze(mean(Md_all,2));
Mi_sub_tt_ss = squeeze(mean(Mi_all,2));
xMX_sub_tt_ss = squeeze(mean(xMX_all,2));
xMd_sub_tt_ss = squeeze(mean(xMd_all,2));
xMi_sub_tt_ss = squeeze(mean(xMi_all,2));

clear MX_all Md_all Mi_all xMX_all xMd_all xMi_all

%%% Mean from tt = ss:end (mean over vertical columns of the sub-diag matrixes)
DT = 0;  %%% I add a DeltaT to avoid the initial spureous (??) peak --------------------> CHECK
MX_sub_ss = zeros(NSUB,Tup);
Md_sub_ss = zeros(NSUB,Tup);
Mi_sub_ss = zeros(NSUB,Tup);
xMX_sub_ss = zeros(NSUB,Tup);
xMd_sub_ss = zeros(NSUB,Tup);
xMi_sub_ss = zeros(NSUB,Tup);
for sub = 1:NSUB
    for ss = 1:Tup-DT
        MX_sub_ss(sub,ss) = TR*trapz(squeeze(MX_sub_tt_ss(sub,ss+DT:end,ss)));
        Md_sub_ss(sub,ss) = TR*trapz(squeeze(Md_sub_tt_ss(sub,ss+DT:end,ss)));
        Mi_sub_ss(sub,ss) = TR*trapz(squeeze(Mi_sub_tt_ss(sub,ss+DT:end,ss)));
        xMX_sub_ss(sub,ss) = TR*trapz(squeeze(xMX_sub_tt_ss(sub,ss+DT:end,ss)));
        xMd_sub_ss(sub,ss) = TR*trapz(squeeze(xMd_sub_tt_ss(sub,ss+DT:end,ss)));
        xMi_sub_ss(sub,ss) = TR*trapz(squeeze(xMi_sub_tt_ss(sub,ss+DT:end,ss)));
    end
end

%%% Mean over different ss values
MX_sub = zeros(1,NSUB);
Md_sub = zeros(1,NSUB);
Mi_sub = zeros(1,NSUB);
xMX_sub = zeros(1,NSUB);
xMd_sub = zeros(1,NSUB);
xMi_sub = zeros(1,NSUB);

for sub = 1:NSUB
    MX_sub(sub) = TR*trapz(squeeze(MX_sub_ss(sub,:)));
    Md_sub(sub) = TR*trapz(squeeze(Md_sub_ss(sub,:)));
    Mi_sub(sub) = TR*trapz(squeeze(Mi_sub_ss(sub,:)));
    xMX_sub(sub) = TR*trapz(squeeze(xMX_sub_ss(sub,:)));
    xMd_sub(sub) = TR*trapz(squeeze(xMd_sub_ss(sub,:)));
    xMi_sub(sub) = TR*trapz(squeeze(xMi_sub_ss(sub,:)));
end

%%% Mean over subjects
MX_final = mean(MX_sub)
Md_final = mean(Md_sub)
Mi_final = mean(Mi_sub)
xMX_final = mean(xMX_sub)
xMd_final = mean(xMd_sub)
xMi_final = mean(xMi_sub)

%%% SAVE FDT METRICS TO FILE %%%
file_metrics = sprintf('metrics_SUB_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d.mat',...
                        CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
save(file_metrics,'MX_sub','Md_sub','Mi_sub','xMX_sub','xMd_sub','xMi_sub', ...
                  'MX_final','Md_final','Mi_final','xMX_final','xMd_final','xMi_final')
%%%%%%%%%%%%%%%%%%%%%%
end %CONDITION
end %NSUBSIM

return

