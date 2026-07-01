function [FC,CV,Cvth,A] = hopf_int(gC,f_diff,sigma)
% HOPF_INT Computes FC and steady-state covariance of the linear Hopf network.
% Author: Juan Manuel Monti
%
% Solves A*Cvth + Cvth*A' = -sigma^2*I (Sylvester equation) for the linearized
% Stuart-Landau network. Bifurcation parameter fixed at a = -0.02.
%
% Inputs:
%   gC     [NxN]  effective connectivity matrix (G*Cij)
%   f_diff [1xN]  per-parcel peak frequencies (Hz)
%   sigma         noise standard deviation
%
% Outputs:
%   FC   [NxN]   functional connectivity (Pearson correlation of x-components)
%   CV   [NxN]   steady-state covariance (x-component block)
%   Cvth [2Nx2N] full steady-state covariance (x and y blocks)
%   A    [2Nx2N] Jacobian of the linearized system
a = -0.02;

N = size(gC,1);
wo = f_diff'*(2*pi);

Cvth = zeros(2*N);

% Jacobian:

s = sum(gC,2);
B = diag(s);

Axx = a*eye(N) - B + gC;
Ayy = Axx;
Axy = -diag(wo);
Ayx = diag(wo);

A = [Axx Axy; Ayx Ayy];
Qn = (sigma^2)*eye(2*N);

Cvth=sylvester(A,A',-Qn);

FCth=corrcov(Cvth);

FC=FCth(1:N,1:N);
CV=Cvth(1:N,1:N);
end
