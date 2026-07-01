function [Ceffsub,FCemp,FCsim] = hopf_linfit_sub(tsdata,C,sub,NPARCELS,TR,f_diff,sigma,Ceffgroup)
% HOPF_LINFIT_SUB Fits subject-specific effective connectivity G*Cij via gradient descent.
% Author: Juan Manuel Monti
%
% Same optimization as hopf_linfit_group but for a single subject, initialized
% from the group-level fit (Ceffgroup). Uses identical hyperparameters:
% Tau=1, epsFC=0.0004, epsFCtau=0.0001, maxC=0.2, up to 5000 iterations.
%
% Inputs:
%   tsdata     [NSUBxNPARCELSxTmax] empirical BOLD time series (all subjects)
%   C          [NxN] structural connectivity (normalized, zero diagonal)
%   sub        subject index (1-based)
%   NPARCELS  number of parcels
%   TR         fMRI repetition time (s)
%   f_diff     [1xN] per-parcel peak frequencies (Hz)
%   sigma      noise standard deviation
%   Ceffgroup  [NxN] group-level effective connectivity (initialization)
%
% Outputs:
%   Ceffsub  [NxN] fitted subject-specific effective connectivity
%   FCemp    [NxN] empirical FC for this subject
%   FCsim    [NxN] simulated FC at convergence

Tau = 1;
epsFC = 0.0004;
epsFCtau = 0.0001;
maxC = 0.2;

SETOLDERROR = 100000;
SETITER = 5000;

indexN = 1:NPARCELS;  %% Cortical areas
N = length(indexN);
% Isubdiag = find(tril(ones(N),-1));

ts = squeeze(tsdata(sub,:,:));

% FC(0)
ts2 = ts(indexN,10:end-10);
% Tm = size(ts2,2);
FCemp = corrcoef(ts2');
% FCPB(sub,:,:) = FCemp;
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
Cnew = Ceffgroup;
olderror = SETOLDERROR;
for iter = 1:SETITER
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
    % errorFC(iter) = mean(mean((FCemp - FCsim).^2));
    % errorCOVtau(iter) = mean(mean((COVtauemp - COVtausim).^2));

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
Ceffsub = Cnew;

% [FCsim,COVsim,COVsimtotal,A] = hopf_int(Ceffsub,f_diff,sigma);
% fittFC_PB(sub) = corr2(FCemp(Isubdiag),FCsim(Isubdiag));
% COVtausim = expm((Tau*TR)*A)*COVsimtotal;
% COVtausim = COVtausim(1:N,1:N);
% for i = 1:N
%     for j = 1:N
%         sigratiosim(i,j) = 1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
%     end
% end
% COVtausim = COVtausim.*sigratiosim;
% fittCVtau_PB(sub) = corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));


% %%% COMPARISON betwwen FCemp(sub) and FCsim(sub)
% fig1 = figure('Visible','off');
% fig1.Position = [100 100 1200 600];
% subplot(1,2,1)
% imagesc(FCemp),colorbar
% title('FCemp sub')
% subplot(1,2,2)
% imagesc(FCsim),colorbar
% title('FCsim sub')
% saveas(fig1,sprintf('2_dset_%d_cond_%d_FCemp_FCsim_sub_%d.png',dataset,condition,sub))
% clear fig1

end
