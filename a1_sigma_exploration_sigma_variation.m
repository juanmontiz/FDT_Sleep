
%%% A1_SIGMA_EXPLORATION_SIGMA_VARIATION  Sweeps sigma to find optimal noise amplitude.
% Author: Juan Manuel Monti
%
% For each CONDITION, runs NSUBSIM nonlinear Hopf simulations per sigma value in
% sigma = 0.01:0.01:0.5 and computes the MSE between empirical and simulated C(t,s)
% (z-score normalized). Requires Ceff_FC_Sleep_COND_*.mat from a0_Ceff_fitting.m.
%
% Error metrics stored as matrices [length(sigma) x NSUBSIM]:
%   errorC_norm_sim_zsc_std   norm(Cemp - Csim)
%   errorC_quad_sim_zsc_std   mean((Cemp - Csim)^2)
%   errorC_d1_sim_zsc_std     mean(|Cemp - Csim|)
%   errorC_d2_sim_zsc_std     sqrt(mean((Cemp - Csim)^2))
%
% Key parameters:
%   NSUBSIM    simulations per sigma value (default: 100)
%   DATAFILTER 0 = unfiltered (default), 1 = filtered
%   HOPFINIT   0 = init once per subject, 1 = init every simulation
%
% Requires: hopf_sim_0init, hopf_sim_1start, make_filt, Ceff_FC_Sleep_COND_*.mat
% Saves:    sigmaexplore_Sleep_COND_<C>_NSUB_<N>_NSUBSIM_<S>_DFILT_<D>_FREQSUB_<F>.mat

for CONDITION = 1:2

clearvars -except CONDITION

%% Choose dataset and other parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATASET
% ./Data/Sleep (15 subjects - 90 parcels)

%%% Number of groups of NSUB
NSUBSIM = 100 %%% <------------ SET THIS

%%% Use filtered or non-filtered data (ts, force, noise) to calculate C, A and R
%   DATAFILTER = 1 --> use filtered data
%   DATAFILTER = 0 --> use non-filtered data
DATAFILTER = 0;

%%% Hopf initialization
%   HOPFINT = 1 --> initializes for each individual simulation
%   HOPFINT = 0 --> initializes ONLY for the first simulation (of each subject)
%                   and all subsequent use the same z0
HOPFINIT = 1;

%%% LINEAR HOPF (not yet implemented)
LIN = 0

%%% Hopf Frequencies
%   FREQSUB = 1 --> different frequencies for each subject
%   FREQSUB = 0 --> frequencies for the mean Power Spectra
FREQSUB = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. LOAD EMPIRIC DATA %%%
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
C=SC/max(max(SC))*0.2;

%%% Parcels
NPARCELS = length(C);
disp(['NPARCELS = ', num2str(NPARCELS)])
Isubdiag = find(tril(ones(NPARCELS),-1));

%%% Subjects in empiric data
NSUB = size(TS_X,2);

disp(['NSUB = ', num2str(NSUB)])
%%% Sets NSUBSIM = NSUB
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

%% Make Filter for signal "MEGstyle" %%%%%%%%%%%
flp = 0.008;            % lowpass frequency of filter
fhi = 0.08;             % highpass
delt = TR;              % sampling interval
k = 2;                  % 2nd order butterworth filter
[bfilt2, afilt2] = make_filt(flp,fhi,delt,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Saves TS data for all subjects up to Tmax %%%
tsdata = zeros(NSUB,NPARCELS,Tmax);
tsdata_filt = zeros(NSUB,NPARCELS,Tmax);
for sub = 1:NSUB
    tsdata(sub,:,:) = TS_X{1,sub}(:,1:Tmax) ;
end
%%% Data filtering
for sub = 1:NSUB
    BOLDdata = (squeeze(tsdata(sub,:,:)));
    timeseriedata = zeros(NPARCELS,Tmax);
    for seed = 1:NPARCELS
        BOLDdata(seed,:) = BOLDdata(seed,:) - mean(BOLDdata(seed,:));
        timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,BOLDdata(seed,:));
    end

    %%% TS empiric filtered and non-filtered (removed mean)
    tsdata_filt(sub,:,:) = timeseriedata;
    tsdata(sub,:,:) = BOLDdata;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  1: NORMALIZATION OF THE TS  %%%
%  N = normalize(___,method,methodtype)
%  "zscore" | "std" (default) | Compute the z-score. Center data to have mean 0, and scale data to have standard deviation 1.
Cemp_sub_zsc_std = zeros(NSUB,NPARCELS,Tmax,Tmax);
%%%
for sub = 1:NSUB
    if DATAFILTER == 0
        tsaux = squeeze(tsdata(sub,:,:));
    elseif DATAFILTER == 1
        tsaux = squeeze(tsdata_filt(sub,:,:));
    end
    for par = 1:NPARCELS
    %%% z-score normalization of the TS
        tsaux_norm(par,:) = normalize(tsaux(par,:),"zscore","std");
        Cemp_sub_zsc_std(sub,par,:,:) = tsaux_norm(par,:)'*tsaux_norm(par,:);
    end
end
clear tsaux tsaux_norm
%%% Group average --> Average over all subjects
Cemp_zsc_std = squeeze(mean(Cemp_sub_zsc_std,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2: Fitting to obtain G.C %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check if fitting was previously done
fileCeff = sprintf('Ceff_FC_Sleep_COND_%d_NSUB_%d_DFILT_%d_FREQSUB_%d.mat',...
                                  CONDITION,NSUB,DATAFILTER,FREQSUB);
if isfile(fileCeff)
    load(fileCeff)
    disp('! Load Ceff from file')
else
    disp('! Run Ceff fitting')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3: Exploration to obtain sigma %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- Exploration to obtain sigma --')
load(fileCeff,'Ceffgroup')

%%% From previous linear fitting
wC = Ceffgroup;
sumC = repmat(sum(wC,2),1,2);

%%% 1 sweep over sig values
sigma = (0.01:0.01:0.5);
for isig = 1:length(sigma)

    %%% MODEL PARAMETERS
    omega1 = repmat(2*pi*f_diff',1,2); omega1(:,1) = -omega1(:,1);
    % Time step
    dt = 0.1*TR/2;
    % noise variance
    sig = sigma(isig)
    
    dsig = sqrt(dt)*sig;
    % "Temperature"
    Temp = sig^2/2;
    
    %%% FIXED ahopf value
    %%% Same as hopf_int.m !!
    ahopf = -0.02;
    
    %%% Frequencies
    omega = omega1;

    %%% 2 Make NSUBSIM "packages" of NSUB
    for sim = 1:NSUBSIM
        sim

        Csub_zsc_std = zeros(NPARCELS,Tup,Tup);
        for sub = 1:NSUB
            % sub
            %%% Initialize Hopf discarding first 2000/dt time steps
            % Initializes EVERY simulation
            if HOPFINIT == 1
                z0 = hopf_sim_0init(ahopf,omega,dsig,dt,NPARCELS,wC,sumC,LIN);
            % Initializes ONLY THE FIRST simulation (all others use the same z0)
            elseif HOPFINIT == 0 && sim == 1
                z0 = hopf_sim_0init(ahopf,omega,dsig,dt,NPARCELS,wC,sumC,LIN);
            end
            %%% Actual Hopf simulation
            [ts, force, noise] = hopf_sim_1start(ahopf,omega,dsig,dt,Tmax,TR,NPARCELS,wC,sumC,LIN,z0);
            %%% Simulated TS
            tserisub = ts;

            %%% Signal filtering %%%
            tssub(:,:) = tserisub(:,Tinf:Tsup);
            tsaux = zeros(NPARCELS,Tmax);
            tsaux_filt = zeros(NPARCELS,Tmax);
            for seed = 1:NPARCELS
                % TS
                tsaux(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:)));
                tsaux_filt(seed,:)  = filtfilt(bfilt2,afilt2,tsaux(seed,:));
                % Hilbert
                Xanalytic = hilbert(demean(tsaux_filt(seed,:)));
                Phases(seed,:) = angle(Xanalytic);
            end
            tssub_filt = tsaux_filt(:,Tinf:Tsup);
            clear ts_aux tsaux_filt

            if DATAFILTER == 0
                tsaux = tssub;
            elseif DATAFILTER == 1
                tsaux = tssub_filt;
            end

            for par = 1:NPARCELS

            %%% z-score normalization of the TS
                tsaux_norm(par,:) = normalize(tsaux(par,:),"zscore","std");
                Csub_zsc_std_aux(par,:,:) = tsaux_norm(par,:)'*tsaux_norm(par,:);
            end
            Csub_zsc_std = (Csub_zsc_std + Csub_zsc_std_aux);
        end
        clear tsaux tsaux_norm
        clear Csub_zsc_std_aux

        %%% Average over NSUB
        Csim_zsc_std = Csub_zsc_std./NSUB;
        clear Csub_zsc_std

        %%% CALCULATE ERRORS %%%%%%%%%%%%%%%%%%%%
        %%% Calculate (empiric - simulated) C(t,s)
        for par = 1:NPARCELS
        %%% z-score Normalization
            % std
            Cempaux = squeeze(Cemp_zsc_std(par,Tinf:Tsup,Tinf:Tsup));
            Csimaux = squeeze(Csim_zsc_std(par,:,:));
            %%% Differences between matrices
            errorC_norm_par_zsc_std(par) = norm(Cempaux - Csimaux);
            errorC_quad_par_zsc_std(par) = mean(mean((Cempaux - Csimaux).^2));
            errorC_d1_par_zsc_std(par) = mean(mean(abs(Cempaux - Csimaux)));
            errorC_d2_par_zsc_std(par) = sqrt(mean(mean((Cempaux - Csimaux).^2)));
        end
        clear Cempaux Csimaux

        %%% Average of the erorr over PARCELS
        %%% z-score Normalization
        errorC_norm_sim_zsc_std(isig,sim) = mean(errorC_norm_par_zsc_std);
        errorC_quad_sim_zsc_std(isig,sim) = mean(errorC_quad_par_zsc_std);
        errorC_d1_sim_zsc_std(isig,sim) = mean(errorC_d1_par_zsc_std);
        errorC_d2_sim_zsc_std(isig,sim) = mean(errorC_d2_par_zsc_std);
    end %SIM
    save(sprintf('sigmaexplore_Sleep_COND_%d_NSUB_%d_NSUBSIM_%d_DFILT_%d_FREQSUB_%d.mat',...
                                     CONDITION,NSUB,NSUBSIM,DATAFILTER,FREQSUB), ...
    'sigma', ...
    'errorC_norm_sim_zsc_std','errorC_quad_sim_zsc_std','errorC_d1_sim_zsc_std','errorC_d2_sim_zsc_std')
end %SIGMA
end %CONDITION

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
