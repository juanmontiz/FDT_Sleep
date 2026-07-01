function [Csimaux,Asimaux,Rsimaux] = funcs_FDT_CAR_sim(tsaux,forceaux,noiseaux,NPARCELS,Tup,Temp)
% FUNCS_FDT_CAR_SIM Computes C(t,s), A(t,s), R(t,s) for a single simulation.
% Author: Juan Manuel Monti
%
% Two-time observables for one Hopf realization (Cugliandolo & Kurchan, 1994).
% Only the lower triangle s <= t is filled (causality constraint).
%
%   C(t,s) = x(t)*x(s)                                 [eq. 2.1]
%   A(t,s) = f(t)*x(s) - f(s)*x(t)                    [eq. 2.10]
%   R(t,s) = (1/2T) * x(t)*eta(s)                     [eq. 2.9]
%
% Inputs:
%   tsaux      [NPARCELSxTup] simulated x-component (BOLD proxy)
%   forceaux   [NPARCELSxTup] deterministic force (negated per eq. 2.7)
%   noiseaux   [NPARCELSxTup] noise divided by sqrt(dt)
%   NPARCELS  number of brain parcels
%   Tup        number of time points (= Tsup - Tinf + 1)
%   Temp       effective temperature (sigma^2 / 2)
%
% Outputs:
%   Csimaux  [NPARCELSxTupxTup] correlation function C(t,s)
%   Asimaux  [NPARCELSxTupxTup] asymmetric correlator A(t,s)
%   Rsimaux  [NPARCELSxTupxTup] linear response function R(t,s)

Csimaux = zeros(NPARCELS,Tup,Tup);
Asimaux = zeros(NPARCELS,Tup,Tup);
Rsimaux = zeros(NPARCELS,Tup,Tup);

%%% First option: nested loops %%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = 1:Tup
   for ss = 1:tt  % must be ss <= tt (causality)
       Csimaux(:,tt,ss) = tsaux(:,tt).*tsaux(:,ss);
       Asimaux(:,tt,ss) = forceaux(:,tt).*tsaux(:,ss) - forceaux(:,ss).*tsaux(:,tt);
       Rsimaux(:,tt,ss) = tsaux(:,tt).*noiseaux(:,ss);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsimaux = (0.5 / Temp) * Rsimaux;
end
