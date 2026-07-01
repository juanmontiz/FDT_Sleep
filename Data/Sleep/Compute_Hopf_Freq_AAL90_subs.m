%%% COMPUTE_HOPF_FREQ_AAL90_SUBS  Computes per-parcel peak frequencies from fMRI power spectra.
% Author: Juan Manuel Monti
%
clear all;

%%% SET CONDITION
% 1 --> Wakefulness (W)
% 2 --> Deep Sleep (N3)
for CONDITION = 1:2;
clearvars -except CONDITION

load('DataSleep_W_N3.mat')
C=SC/max(max(SC))*0.2;
NPARCELS = length(C);

%%% TS
if CONDITION == 1
    TS_X = TS_W;
elseif CONDITION == 2
    TS_X = TS_N3;
end

%%% Subjects and Tmax in empiric data
NSUB = size(TS_X,2);
for i = 1:size(TS_X,2); Tmax(i) = size(TS_X{1,i},2); end
Tmax = min(min(Tmax));

%%%
% Parameters of the data
TR=2;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELS),-1));

%%%%
fce1=zeros(NSUB,NPARCELS,NPARCELS);

for sub=1:NSUB
    sub

    ts = TS_X{1,sub}(:,1:Tmax);
    
    % [Ns, Tmax]=size(ts);
    TT=Tmax;
    tss=zeros(NPARCELS,Tmax);
    Ts = TT*TR;
    % Possible frequencies for the nodes
    freq = (0:TT/2-1)/Ts;
    nfreqs = length(freq);
    
    % Process the time series: detrend-mean ; filtfilt ; fft
    for seed = 1:NPARCELS
        x(seed,:) = detrend(ts(seed,:) - mean(ts(seed,:)));
        tss(seed,:) = filtfilt(bfilt,afilt,x(seed,:));
        pw = abs(fft(tss(seed,:)));
        % Power spectrum
        PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    fce1(sub,:,:) = corrcoef(tss','rows','pairwise');
end
fce_sub = fce1;
% Mean functional connectivity empirical (averages fce1 over all subjects)
fce=squeeze(mean(fce1));

%%%%% Frequencies for each subject %%%
for sub = 1:NSUB
    Power_Areas=squeeze(PowSpect(:,:,sub));
    for seed = 1:NPARCELS
        Power_Areas_gfilt(:,seed) = gaussfilt(freq,Power_Areas(:,seed)',0.01);
    end

    % Sets the freq. of each parcel by searching those that give a maximum
    % in the power spectra for each parcel
    [maxpowdata,index]=max(Power_Areas_gfilt);
    f_diff = freq(index);
    % If any value of freq. is == 0 then takes the average of all the
    % non-zero values
    f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

    f_diff_sub(sub,:) = f_diff;
end

%%%%% Frequencies for the mean power spectra %%%
%     mean(PowSpect,3)-> averages PowSpect over all subjects
%     then applies Gauss filter
Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:NPARCELS
    Power_Areas_gfilt(:,seed) = gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

% Sets the freq. of each parcel by searching those that give a maximum
% in the power spectra for each parcel
[maxpowdata,index]=max(Power_Areas_gfilt);
f_diff = freq(index);
% If any value of freq. is == 0 then takes the average of all the
% non-zero values
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

if CONDITION == 1
    file_save = 'hopf_freq_AAL90_COND_1_W.mat';
elseif CONDITION == 2
    file_save = 'hopf_freq_AAL90_COND_2_N3.mat';
end
save(file_save,'f_diff','f_diff_sub','fce','fce_sub')

end
