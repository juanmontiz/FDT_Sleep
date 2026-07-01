
%%% A0_CEFF_FITTING  Fits effective connectivity (G*Cij) from empirical fMRI data.
% Author: Juan Manuel Monti
%
% For each CONDITION (1=Wakefulness, 2=Deep Sleep N3):
%   1. Loads SC and BOLD time series from Data/Sleep/DataSleep_W_N3.mat
%   2. Loads per-parcel oscillation frequencies from hopf_freq_AAL90_COND_*.mat
%   3. Fits group-level Ceff via gradient descent on FC(0) and COV(tau=1)
%   4. Fits per-subject Ceff initialized from group fit (parfor over subjects)
%   5. Saves to Ceff_FC_Sleep_COND_<C>_NSUB_<N>_DFILT_<D>_FREQSUB_<F>.mat
%
% Key parameters:
%   DATAFILTER  0 = unfiltered (default), 1 = bandpass filtered
%   FREQSUB     0 = group-mean frequencies (default), 1 = subject-specific
%   SETSIGMA    noise std used in hopf_int: 0.12 (W), 0.06 (N3)
%
% Requires: hopf_linfit_group, hopf_linfit_sub, hopf_int, make_filt
% Saves:    Ceff_FC_Sleep_COND_*.mat

conds = [1 2];

%%% BRAIN STATES
% CONDITION = 1 --> Wakefulness (W)
% CONDITION = 2 --> Deep Sleep (N3)
for CONDITION = conds
CONDITION

clearvars -except DATASET CONDITION FREQSUB

%% Choose dataset and other parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATASET
% ./Data/Sleep [Wakefulness (15 subs.) ; Sleep (15 subs.) - 90 AAL parcels]

%%% Use filtered or non-filtered data (ts, force, noise) to calculate C, A and R
% DATAFILTER = 1 --> use filtered data
% DATAFILTER = 0 --> use non-filtered data
DATAFILTER = 0;

%%% Set value of standard deviation of noise
if CONDITION == 1
    SETSIGMA = 0.12;
elseif CONDITION == 2
    SETSIGMA = 0.06;
end

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


%% 1. LOAD EMPIRIC DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- 1. Load empirical data --')
% Structural Connectivity and Time series
load('./Data/Sleep/DataSleep_W_N3.mat')
% Nodes Frequencies and FC empiric, the output of Compute_Hopf_freq_AAL90
if CONDITION == 1
    load('./Data/Sleep/hopf_freq_AAL90_COND_1_W.mat') % the output of Compute_Hopf_freq_AAL90
    %%% Time-series (TS)
    TS_X = TS_W;
elseif  CONDITION == 2
    load('./Data/Sleep/hopf_freq_AAL90_COND_2_N3.mat') % the output of Compute_Hopf_freq_AAL90
    %%% Time-series
    TS_X = TS_N3;
end
%%% Structural Connectivity (C)
C=SC/max(max(SC))*0.2;
% Diagonal of C equal zero
C = C.*~eye(size(C));

%%% Parcels
NPARCELS = length(C);
disp(['NPARCELS = ', num2str(NPARCELS)])

%%% Subjects in empiric data
NSUB = size(TS_X,2);
disp(['NSUB = ', num2str(NSUB)])

%%% Tmax & TR
Tmax = zeros(1,NSUB);
for i = 1:NSUB; Tmax(i) = size(TS_X{1,i},2); end
Tmax = min(min(Tmax));
disp(['Tmax = ', num2str(Tmax)])
Tinf = 10;
Tsup = Tmax - Tinf;
Tup = Tsup - Tinf + 1;
disp(['Tinf = ', num2str(Tinf)])
disp(['Tsup = ', num2str(Tsup)])
% Repetition time
TR = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Make Filter for signal "MEGstyle" %%%%%%%%%%%
flp = 0.008;            % lowpass frequency of filter
fhi = 0.08;             % highpass
delt = TR;              % sampling interval
k = 2;                  % 2nd order butterworth filter
[bfilt2, afilt2] = make_filt(flp,fhi,delt,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TS data for all subjects up to Tmax %%%%%%%%%
tsdata = zeros(NSUB,NPARCELS,Tmax);
tsdata_filt = zeros(NSUB,NPARCELS,Tmax);
for sub = 1:NSUB
    tsdata(sub,:,:) = TS_X{1,sub}(:,1:Tmax) ;
end
%%% Data filtering
for sub = 1:NSUB
    BOLDdata = (squeeze(tsdata(sub,:,:)));
    timeseriedata = zeros(NPARCELS,Tmax);
    for np = 1:NPARCELS
        BOLDdata(np,:) = BOLDdata(np,:) - mean(BOLDdata(np,:));
        timeseriedata(np,:) = filtfilt(bfilt2,afilt2,BOLDdata(np,:));
    end

    %%% TS empiric filtered and non-filtered (removed mean)
    tsdata_filt(sub,:,:) = timeseriedata;
    tsdata(sub,:,:) = BOLDdata;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. Fitting to obtain Ceff (Effective Connectivity) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- 2. Fitting to obtain Ceff = G.Cij --')
%%% Check if fitting was previously done
fileCeff = sprintf('Ceff_FC_Sleep_COND_%d_NSUB_%d_DFILT_%d_FREQSUB_%d.mat',...
                                  CONDITION,NSUB,DATAFILTER,FREQSUB);
if isfile(fileCeff)
    load(fileCeff)
    dolinfit = 0;
    disp('! Load Ceff from file')
else
    dolinfit = 1;
end
sigma = 0.05; % Not at all sensitive to sigma but needs a value to work %
%%% Fitting of EC using a Linear Hopf
if dolinfit == 1
    if DATAFILTER == 0
        tsaux = tsdata;
    elseif DATAFILTER == 1
        tsaux = tsdata_filt;
    end
    disp('! Do Fitting to obtain Ceff')
    %%% Group EC
    disp('+: Linear fit for GROUP')
    [Ceffgroup,FCempgroup,FCsimgroup] = hopf_linfit_group(tsaux,C,NSUB,NPARCELS,TR,f_diff,sigma);
    %%% Individual EC for each subject
    Ceffsub = zeros(NSUB,NPARCELS,NPARCELS);
    FCempsub = zeros(NSUB,NPARCELS,NPARCELS);
    FCsimsub = zeros(NSUB,NPARCELS,NPARCELS);
    parfor sub = 1:NSUB
        disp(['+: Linear fit for sub. = ', num2str(sub), ' of ', num2str(NSUB)])
        if FREQSUB == 0
            f_diff_aux = f_diff;
        elseif FREQSUB == 1
            f_diff_aux = squeeze(f_diff_sub(sub,:));
        end
        [Ceffsub(sub,:,:),FCempsub(sub,:,:),FCsimsub(sub,:,:)] = hopf_linfit_sub(tsaux,C,sub,NPARCELS,TR,f_diff_aux,sigma,Ceffgroup);
    end
    clear f_diff_aux
    %%% Save the empiric EC and FC to use later
    save(fileCeff,'Ceffsub','Ceffgroup','FCempgroup','FCsimgroup','FCempsub','FCsimsub')
end
clear tsaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %CONDITION
return
