function [Ceffgroup,FCemp,FCsim] = hopf_linfit_group(tsdata,C,NSUB,NPARCELS,TR,f_diff,sigma)
% HOPF_LINFIT_GROUP Fits the group-averaged effective connectivity G*Cij via gradient descent.
% Author: Juan Manuel Monti
%
% Minimizes combined MSE of FC(0) and COV(tau=1) between empirical (group-averaged)
% and linear Hopf model. Cij updated only where SC(i,j)>0 or j is the contralateral
% mirror (j == N-i+1). Fixed hyperparameters: Tau=1, epsFC=0.0004, epsFCtau=0.0001,
% maxC=0.2, maxIter=5000. Stops early if relative improvement < 0.1% or error rises.
%
% Inputs:
%   tsdata     [NSUBxNPARCELSxTmax] empirical BOLD time series
%   C          [NxN] structural connectivity (normalized, zero diagonal)
%   NSUB       number of subjects
%   NPARCELS  number of parcels
%   TR         fMRI repetition time (s)
%   f_diff     [1xN] per-parcel peak frequencies (Hz)
%   sigma      noise standard deviation (passed to hopf_int)
%
% Outputs:
%   Ceffgroup  [NxN] fitted group-level effective connectivity
%   FCemp      [NxN] empirical FC averaged over subjects
%   FCsim      [NxN] simulated FC at convergence

Tau = 1;
epsFC=0.0004;
epsFCtau=0.0001;
maxC=0.2;

indexN = 1:NPARCELS;  %% Cortical areas
N = length(indexN);

for sub = 1:NSUB
    ts = squeeze(tsdata(sub,:,:));
    % FC(0)
    ts2 = ts(indexN,10:end-10);
    Tm = size(ts2,2);
    FCemp = corrcoef(ts2');
    FCPB(sub,:,:) = FCemp;
    COVemp = cov(ts2');
    % COV(tau)
    tst = ts2';
    for i = 1:N
        for j = 1:N
            sigratio(i,j) = 1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag, lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx = find(lags == Tau);
            COVtauemp(i,j) = clag(indx)/size(tst,1);
        end
    end
    COVtauemp = COVtauemp.*sigratio;
    COVtauPB(sub,:,:) = COVtauemp;
end
FCemp = squeeze(mean(FCPB));
COVtauemp = squeeze(mean(COVtauPB));
    
    
% %%% COMPARISON betwwen FCemp0 and FCemp (mean already taken)
% fig1 = figure('Visible','off');
% fig1.Position = [100 100 1200 600];
% subplot(1,2,1)
% imagesc(FCemp0),colorbar
% title('FCemp 0')
% subplot(1,2,2)
% imagesc(FCemp),colorbar
% title('FCemp linfit')
% saveas(fig1,sprintf('1_dset_%d_cond_%d_FCemp0_FCemp.png',dataset,condition))
% clear fig1

Cnew = C;
olderror = 100000;
for iter = 1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A] = hopf_int(Cnew,f_diff,sigma);
    COVtausim = expm((Tau*TR)*A)*COVsimtotal;
    COVtausim = COVtausim(1:N,1:N);
    for i = 1:N
        for j = 1:N
            sigratiosim(i,j) = 1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim = COVtausim.*sigratiosim;
    errorFC(iter) = mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter) = mean(mean((COVtauemp-COVtausim).^2));

    if mod(iter,100) < 0.1
        errornow = mean(mean((FCemp - FCsim).^2)) + mean(mean((COVtauemp - COVtausim).^2));
        if  (olderror - errornow)/errornow < 0.001
            break;
        end
        if  olderror < errornow
            break;
        end
        olderror = errornow;
    end

    for i = 1:N  %% learning
        for j = 1:N
            if (C(i,j) > 0 || j == N-i+1)
                Cnew(i,j) = Cnew(i,j) + epsFC*(FCemp(i,j) - FCsim(i,j)) ...
                                      + epsFCtau*(COVtauemp(i,j) - COVtausim(i,j));
                if Cnew(i,j) < 0
                    Cnew(i,j) = 0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end    
Ceffgroup = Cnew;
end

