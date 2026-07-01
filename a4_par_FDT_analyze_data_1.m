%%% A4_PAR_FDT_ANALYZE_DATA_1  Reduces FDT matrices to per-parcel scalar metrics.
% Author: Juan Manuel Monti
%
% Loads partial_Sleep_COND_*.mat from a3 and integrates FDT observables over (t,s)
% via trapz to obtain one scalar per parcel (averaged over subjects first).
%
% Reduction steps:
%   1. Average |XFDR-1|, |dVFDT|, |iVFDT| over subjects  -> [NPARCELSxTupxTup]
%   2. Integrate over t (trapz) for each fixed s           -> [NPARCELSxTup]
%   3. Integrate over s (trapz)                            -> [NPARCELSx1]
%
% Output metrics (one scalar per parcel):
%   MX_par   |XFDR - 1| integrated over (t,s)
%   Md_par   |dVFDT|    integrated over (t,s)
%   Mi_par   |iVFDT|    integrated over (t,s)
%   xMX_par, xMd_par, xMi_par  same using alternative response (xR)
%
% Requires: partial_Sleep_COND_*.mat
% Saves:    metrics_PAR_Sleep_COND_<C>_NSUB_<N>_NSIM_<S>_DFILT_<D>_SIGMA_<SG>_LIN_<L>_HINIT_<H>_FREQSUB_<F>_trapz.mat

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

%%% Mean over SUBS
MX_par_tt_ss = squeeze(mean(MX_all,1));
Md_par_tt_ss = squeeze(mean(Md_all,1));
Mi_par_tt_ss = squeeze(mean(Mi_all,1));
xMX_par_tt_ss = squeeze(mean(xMX_all,1));
xMd_par_tt_ss = squeeze(mean(xMd_all,1));
xMi_par_tt_ss = squeeze(mean(xMi_all,1));

clear MX_all Md_all Mi_all xMX_all xMd_all xMi_all

%%% Integral from tt = ss:end (integration over vertical columns of the sub-diag matrixes)
DT = 0;  %%% I add a DeltaT to avoid the initial spureous peak (?) --------------------> CHECK
MX_par_ss = zeros(NPARCELS,Tup);
Md_par_ss = zeros(NPARCELS,Tup);
Mi_par_ss = zeros(NPARCELS,Tup);
xMX_par_ss = zeros(NPARCELS,Tup);
xMd_par_ss = zeros(NPARCELS,Tup);
xMi_par_ss = zeros(NPARCELS,Tup);
for par = 1:NPARCELS
    for ss = 1:Tup-DT
        MX_par_ss(par,ss) = TR*trapz(squeeze(MX_par_tt_ss(par,ss+DT:end,ss)));
        Md_par_ss(par,ss) = TR*trapz(squeeze(Md_par_tt_ss(par,ss+DT:end,ss)));
        Mi_par_ss(par,ss) = TR*trapz(squeeze(Mi_par_tt_ss(par,ss+DT:end,ss)));
        xMX_par_ss(par,ss) = TR*trapz(squeeze(xMX_par_tt_ss(par,ss+DT:end,ss)));
        xMd_par_ss(par,ss) = TR*trapz(squeeze(xMd_par_tt_ss(par,ss+DT:end,ss)));
        xMi_par_ss(par,ss) = TR*trapz(squeeze(xMi_par_tt_ss(par,ss+DT:end,ss)));
    end
end

%%% Mean over different ss values
MX_par = zeros(1,NPARCELS);
Md_par = zeros(1,NPARCELS);
Mi_par = zeros(1,NPARCELS);
xMX_par = zeros(1,NPARCELS);
xMd_par = zeros(1,NPARCELS);
xMi_par = zeros(1,NPARCELS);

for par = 1:NPARCELS
    MX_par(par) = TR*trapz(squeeze(MX_par_ss(par,:)));
    Md_par(par) = TR*trapz(squeeze(Md_par_ss(par,:)));
    Mi_par(par) = TR*trapz(squeeze(Mi_par_ss(par,:)));
    xMX_par(par) = TR*trapz(squeeze(xMX_par_ss(par,:)));
    xMd_par(par) = TR*trapz(squeeze(xMd_par_ss(par,:)));
    xMi_par(par) = TR*trapz(squeeze(xMi_par_ss(par,:)));
end

%%% SAVE FDT METRICS TO FILE %%%
file_metrics = sprintf('metrics_PAR_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d_trapz.mat',...
                        CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);
save(file_metrics,'MX_par','Md_par','Mi_par','xMX_par','xMd_par','xMi_par')
%%%%%%%%%%%%%%%%%%%%%%
end %CONDITION
end %NSUBSIM

return

%  imagesc(squeeze(Mi_sub_tt_ss(1,:,:))), ylabel('t'), xlabel('s'), colorbar
%  plot(squeeze(Mi_sub_tt_ss(1,:,1))), ylabel('Mi(1,:,1)'), xlabel('s')
%  plot(squeeze(Mi_sub_tt_ss(1,:,35))), ylabel('Mi(1,:,35)'), xlabel('s')
%  plot(squeeze(Mi_sub_ss(1,:))), ylabel('Mi(1,:)'), xlabel('s')
