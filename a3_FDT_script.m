%%% A3_FDT_SCRIPT  Computes C(t,s), A(t,s), R(t,s) and FDT violation metrics via Hopf simulation.
% Author: Juan Manuel Monti
%
% For each CONDITION and NSUBSIM (convergence sweep [100 250 750 1000 2500 7500 10000]),
% runs NSUBSIM nonlinear Hopf simulations per subject and computes two-time observables
% following Cugliandolo & Kurchan (1994). Saves partial results after each subject.
%
% Computed quantities (per subject x parcel):
%   Csub     correlation function  C(t,s) = <x(t)*x(s)>          [eq. 2.1]
%   Asub     asymmetric correlator A(t,s)                         [eq. 2.10]
%   Rsub     response function     R(t,s) = <x(t)*eta(s)>/(2T)   [eq. 2.9]
%   XFDRsub  FDT ratio X = T*R/(dC/ds)                           [Lippiello 1999 eq. 43]
%   dVFDTsub differential violation dC/ds - T*R                   [Cugliandolo 1997 eq. 1]
%   iVFDTsub integral violation C(t,t)-C(t,s)-T*int(R,ss:tt)     [Cugliandolo 1997 eq. 2]
%   xR*, xXFDR*, xdVFDT*, xiVFDT*: alternative versions using R from C derivatives and A
%
% Key parameters:
%   SETSIGMA   0.12 (W), 0.06 (N3)
%   HOPFINIT   0 = init once per subject, 1 = init every simulation
%   DATAFILTER 0 = unfiltered (default)
%
% Requires: hopf_sim_0init, hopf_sim_1start, funcs_FDT_CAR_sim, derivative
%           Ceff_FC_Sleep_COND_*.mat
% Saves:    partial_Sleep_COND_<C>_NSUB_<N>_NSIM_<S>_DFILT_<D>_SIGMA_<SG>_LIN_<L>_HINIT_<H>_FREQSUB_<F>.mat

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

%%% Number of simulations for each subject
Nsimulations = [100 250 750 1000 2500 7500 10000]; % --> Use different number of simulations to check convergence of metrics
% Nsimulations = 10000; % --> 10k simulations is OK
for NSUBSIM = Nsimulations
NSUBSIM

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

%%% Linearization:
%   LIN = 0 --> Cubic term in Hopf model
%   LIN = 1 --> Linear Hopf model (NOT YET IMPLEMENTED.!)
LIN = 0;

%%% Hopf initialization
%   HOPFINT = 1 --> initializes for each individual simulation
%   HOPFINT = 0 --> initializes ONLY for the first simulation (of each subject)
%                   and all subsequent use the same z0
HOPFINIT = 1;

%%% Hopf Frequencies
%   FREQSUB = 1 --> different frequencies for each subject
%   FREQSUB = 0 --> frequencies for the mean Power Spectra
FREQSUB = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. LOAD EMPIRIC DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- 1. Load empirical data --')
% Structural Connectivity and Time series
%  load('./Data/Sleep/DataSleep_W_N3.mat')
% Nodes Frequencies and FC empiric, the output of Compute_Hopf_freq_AAL90
if CONDITION == 1
    load('./Data/Sleep/hopf_freq_AAL90_COND_1_W.mat') % the output of Compute_Hopf_freq_AAL90
    %%% TS
    TS_X = TS_W;
elseif  CONDITION == 2
    load('./Data/Sleep/hopf_freq_AAL90_COND_2_N3.mat') % the output of Compute_Hopf_freq_AAL90
    %%% TS
    TS_X = TS_N3;
end
%%% SC
%  C=SC/max(max(SC))*0.2;
%%% Diagonal of C equal zero
%  C = C.*~eye(size(C));

%%% Parcels
%  NPARCELS = length(C);
NPARCELS = 90;
disp(['NPARCELS = ', num2str(NPARCELS)])

%%% Subjects in empiric data
NSUB = size(TS_X,2);
disp(['NSUB = ', num2str(NSUB)])

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

%%% Make Filter for signal "MEGstyle" %%%%%%%%%%%
flp = 0.008;            % lowpass frequency of filter
fhi = 0.08;             % highpass
delt = TR;              % sampling interval
k = 2;                  % 2nd order butterworth filter
[bfilt2, afilt2] = make_filt(flp,fhi,delt,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Saves TS data for all subjects up to Tmax %%%
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


%% 2: Load Effective Connectivity %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- 2. Effective Conectivity --')
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


%% 3: Hopf Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- 3. Hopf Model --')

%%% MODEL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step
dt = 0.1*TR/2;
% noise variance
sig = SETSIGMA;
dsig = sqrt(dt)*sig;
% "Temperature"
Temp = sig^2/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIXED ahopf value %%%
%%% Same as hopf_int.m !! %%%
ahopf = -0.02;
disp(['Calculation for fixed ahopf = ', num2str(ahopf)])

%%% Coupling and Connectivity %%%
%%% This is fitted earlier (G*C --> Effective Connectivity)
% wC = G*C;

%%% Frequencies %%%
omega1 = repmat(2*pi*f_diff',1,2); omega1(:,1) = -omega1(:,1);
omega = omega1;

%%% Initialize subject matrices
FCsim_sub_filt = zeros(NSUB,NPARCELS,NPARCELS);
FCsim_sub = zeros(NSUB,NPARCELS,NPARCELS);
Csub = zeros(NSUB,NPARCELS,Tup,Tup);
Asub = zeros(NSUB,NPARCELS,Tup,Tup);
Rsub = zeros(NSUB,NPARCELS,Tup,Tup);
xRsub = zeros(NSUB,NPARCELS,Tup,Tup);
dtCsub = zeros(NSUB,NPARCELS,Tup,Tup);
dsCsub = zeros(NSUB,NPARCELS,Tup,Tup);
iVFDTsub = zeros(NSUB,NPARCELS,Tup,Tup);
xiVFDTsub = zeros(NSUB,NPARCELS,Tup,Tup);

%%% Check for previous partial calculation
substart = 1;
file_partial = sprintf('partial_Sleep_COND_%d_NSUB_%d_NSIM_%d_DFILT_%d_SIGMA_%0.2f_LIN_%d_HINIT_%d_FREQSUB_%d.mat',...
                                      CONDITION,NSUB,NSUBSIM,DATAFILTER,SETSIGMA,LIN,HOPFINIT,FREQSUB);

%%% Model for each subject %%%
%%% Start simulations for each subject
for sub = substart:NSUB
    disp(['+: Subject = ', num2str(sub), ' of ', num2str(NSUB)])

    %%% Option for considering frequencies for each sub
    if FREQSUB == 1
        f_diff_aux = squeeze(f_diff_sub(sub,:));
        omega1 = repmat(2*pi*f_diff_aux',1,2); omega1(:,1) = -omega1(:,1);
        omega = omega1;
    end

    %%% wC matrix (G.Cij) is obtained by fitting the FC and COVtau
    wC = squeeze(Ceffsub(sub,:,:));
    sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj

    %%% Initialize matrices used in simulations
    tssim = zeros(NPARCELS,Tup);
    fosim = zeros(NPARCELS,Tup);
    nosim = zeros(NPARCELS,Tup);
    FCsim_filt = zeros(NPARCELS,NPARCELS);
    FCsim = zeros(NPARCELS,NPARCELS);
    Csim = zeros(NPARCELS,Tup,Tup);
    Asim = zeros(NPARCELS,Tup,Tup);
    Rsim = zeros(NPARCELS,Tup,Tup);
    %%% I run NSUBSIM for each subject
    for sim = 1:NSUBSIM

        %%% Initialize Hopf discarding first 2000/dt time steps:
        % Initializes EVERY simulation || Initializes ONLY THE FIRST simulation (all others use the same z0)
        if HOPFINIT == 1 || (HOPFINIT == 0 && sim == 1)
            z0 = hopf_sim_0init(ahopf,omega,dsig,dt,NPARCELS,wC,sumC,LIN);
        end
        %%% Actual Hopf simulation        
        [ts, force, noise] = hopf_sim_1start(ahopf,omega,dsig,dt,Tmax,TR,NPARCELS,wC,sumC,LIN,z0);

        %%% Simulated TS, FORCE and NOISE for each simulation
        %%% Signal filtering %%%
        ts_aux = zeros(NPARCELS,Tmax);
        ts_filt_aux = zeros(NPARCELS,Tmax);
        fo_aux = zeros(NPARCELS,Tmax);
        fo_filt_aux = zeros(NPARCELS,Tmax);
        no_aux = zeros(NPARCELS,Tmax);
        no_filt_aux = zeros(NPARCELS,Tmax);
        for np = 1:NPARCELS
            % TS
            ts_aux(np,:) = detrend(ts(np,:) - mean(ts(np,:)));
            ts_filt_aux(np,:)  = filtfilt(bfilt2,afilt2,ts_aux(np,:));
            % Force
            fo_aux(np,:) = detrend(force(np,:) - mean(force(np,:)));
            fo_filt_aux(np,:)  = filtfilt(bfilt2,afilt2,fo_aux(np,:));
            % Noise
            no_aux(np,:) = detrend(noise(np,:) - mean(noise(np,:)));
            no_filt_aux(np,:)  = filtfilt(bfilt2,afilt2,no_aux(np,:));
        end
        clear ts force noise
        %
        tssim = ts_aux(:,Tinf:Tsup);
        fosim = fo_aux(:,Tinf:Tsup);
        nosim = no_aux(:,Tinf:Tsup);
        tssim_filt = ts_filt_aux(:,Tinf:Tsup);
        fosim_filt = fo_filt_aux(:,Tinf:Tsup);
        nosim_filt = no_filt_aux(:,Tinf:Tsup);
        clear ts_aux fo_aux no_aux ts_filt_aux fo_filt_aux no_filt_aux
        %
        if DATAFILTER == 0
            tsaux = tssim;
            forceaux = fosim(:,:);
            noiseaux = nosim(:,:);
        elseif DATAFILTER == 1
            tsaux = tssim_filt;
            forceaux = fosim_filt(:,:);
            noiseaux = nosim_filt(:,:);
        end

        %%% FC Simulated
        FCsim = FCsim + corrcoef(tssim');
        FCsim_filt = FCsim_filt + corrcoef(tssim_filt');

        %%% Calculates C(t,s), A(t,s) and R(t,s) for one simulation
        noiseaux = noiseaux/sqrt(dt); % We must divide by sqrt(dt)
        [Csimaux,Asimaux,Rsimaux] = funcs_FDT_CAR_sim(tsaux,forceaux,noiseaux,NPARCELS,Tup,Temp);

        %%% Sums over all simulations
        Csim = (Csim + Csimaux);
        Asim = (Asim + Asimaux);
        Rsim = (Rsim + Rsimaux);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %sim
    clear tsaux forceaux noiseaux Csimaux Asimaux Rsimaux

    %%% AVERAGING FC OVER SIMULATIONS
    FCsim_sub(sub,:,:) = FCsim./NSUBSIM;
    FCsim_sub_filt(sub,:,:) = FCsim_filt./NSUBSIM;
    clear FCsim_filt FCsim

    %%% <EMSEMBLE AVERAGE>
    Csub(sub,:,:,:) = Csim./NSUBSIM;
    Asub(sub,:,:,:) = Asim./NSUBSIM;
    Rsub(sub,:,:,:) = Rsim./NSUBSIM;
    clear Csim Asim Rsim

    %%% Differentiation of C in t and s
    dtCaux = zeros(NPARCELS,Tup,Tup);
    dsCaux = zeros(NPARCELS,Tup,Tup);
    dxaux = TR*(1:Tup);
    for np = 1:NPARCELS
        Caux = squeeze(Csub(sub,np,:,:));
        for ti = 1:Tup
            dtCaux(np,:,ti) = derivative(dxaux,Caux(:,ti));
            dsCaux(np,ti,:) = derivative(dxaux,Caux(ti,:));
            %%% Uncomment to check differentiation: %%%
            % Ctest_t(np,:,ti) = cumtrapz(dxaux,squeeze(dtCaux(np,:,ti))) + Caux(1,ti);
            % Ctest_s(np,ti,:) = cumtrapz(dxaux,squeeze(dsCaux(np,ti,:))) + Caux(ti,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    dtCsub(sub,:,:,:) = dtCaux;
    dsCsub(sub,:,:,:) = dsCaux;
    %%% Uncomment to check differentiation: %%%
    % Ctest_t_sub(sub,:,:,:) = Ctest_t;
    % Ctest_s_sub(sub,:,:,:) = Ctest_s;
    % plot(squeeze(Csub(1,1,:,1)),'LineWidth',1.5), hold on,
    % plot(squeeze(Ctest_t_sub(1,1,:,1)),'LineWidth',1.5),
    % plot(squeeze(Ctest_s_sub(1,1,:,1)),'--','LineWidth',1.5), hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear Caux dtCaux dsCaux

    %%% R(t,s) as in Cugliandolo 1994 eq. (2.9)
    xRsub = (1.0/2.0/Temp) * (dsCsub - dtCsub - Asub);

    %%% FDR as in Lippiello 1999 eq. (43)
    xXFDRsub = 0.5 * (1 - dtCsub ./ dsCsub) - 0.5 * Asub ./ dsCsub;
    % I remove NANs that may appear when dividing by zero-values of dsCsim
    xXFDRsub(isnan(xXFDRsub)) = 0;
    %%% OBS: also could be calculated as:
    %   X = Temp * xRsub / dsCsub
    %     = Temp * Rsub / dsCsub
    XFDRsub = Temp * Rsub ./ dsCsub;

    %%% dVFDT Differential Violation of FDT
    %   as defined in Cugliandolo 1997 eq. (1)
    %   dVFDT = dC(t,s)/ds - T * R(t,s)
    dVFDTsub = dsCsub - Temp * Rsub;
    xdVFDTsub = dsCsub - Temp * xRsub;

    %%% iVFDT Integral violation of FDT
    %%% as defined in Cugliandolo 1997 eq. (2)
    iVFDTsub(sub,:,1,1) = 0;
    xiVFDTsub(sub,:,1,1) = 0;
    iVFDTsub(sub,:,Tup,Tup) = 0;
    xiVFDTsub(sub,:,Tup,Tup) = 0;
    for np = 1:NPARCELS
        for tt = 2:Tup
            for ss = 1:tt-1
                tintaux = TR * (ss:tt);
                %
                Rintaux = squeeze(Rsub(sub,np,tt,(ss:tt)));
                intRaux = trapz(tintaux,Rintaux);
                iVFDTsub(sub,np,tt,ss) = Csub(sub,np,tt,tt) - Csub(sub,np,tt,ss) - Temp * intRaux;
                %
                Rintaux = squeeze(xRsub(sub,np,tt,(ss:tt)));
                intRaux = trapz(tintaux,Rintaux);
                xiVFDTsub(sub,np,tt,ss) = Csub(sub,np,tt,tt) - Csub(sub,np,tt,ss) - Temp * intRaux;
            end
        end
    end
    clear tintaux intRaux

    %%% HERE I SAVE THE PARTIAL PROGRESS UP TO THIS SUBJECT TO CONTINUE AFTER
    sublast = sub;
    saveOK = 0;
    save(file_partial,'saveOK')
    save(file_partial,'sublast','-append')
    save(file_partial,'saveOK','FCsim_sub','FCsim_sub_filt',...
                      'Csub','Asub','Rsub','xRsub','dtCsub','dsCsub',...
                      'XFDRsub','xXFDRsub','dVFDTsub','xdVFDTsub','iVFDTsub','xiVFDTsub',...
                      'Temp','sublast') 
    saveOK = 1;
    save(file_partial,'saveOK','-append')
end %sub
end %NSUBSIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %CONDITION
return
