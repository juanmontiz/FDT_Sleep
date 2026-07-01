function z0 = hopf_sim_0init(ahopf,omega,dsig,dt,NPARCELS,wC,sumC,LIN)
% HOPF_SIM_0INIT Burn-in phase of Hopf simulation to obtain an equilibrated initial state.
% Author: Juan Manuel Monti
%
% Integrates the Stuart-Landau network with Euler-Maruyama for 2000 s (2000/dt steps),
% returning the final state z0 for use as the starting point of recorded simulations.
%
% Inputs:
%   ahopf      bifurcation parameter (fixed: -0.02)
%   omega      [NPARCELSx2] angular frequencies (rad/s); column 1 negated
%   dsig       noise amplitude pre-scaled by sqrt(dt)
%   dt         integration timestep (s)
%   NPARCELS  number of brain parcels
%   wC         [NPARCELSxNPARCELS] effective connectivity (G*Cij)
%   sumC       [NPARCELSx2] row sums of wC replicated for x and y components
%   LIN        0 = cubic (nonlinear) Hopf; 1 = linear (not yet implemented)
%
% Output:
%   z0  [NPARCELSx2] equilibrated state [x, y] after burn-in

%%% Hopf Simulation
a = ahopf.*ones(NPARCELS,2);
z = 0.1*ones(NPARCELS,2); % Initialize z --> x = z(:,1), y = z(:,2)

%%% Discard first 2000/dt time steps
for t = 0:dt:2000
    %%% Coupling and "flipped" z
    sum = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    %%%%%% Hopf Terms %%%%%%
    %%% Noise term in Hopf
    % Different noise value for Re and Im
    no = dsig*randn(NPARCELS,2);
    % % % Same noise value for Re and Im
    % % noise = dsig*randn(NPARCELS,1);
    % % no = [noise,noise];
    %%% Force term in Hopf
    % Cubic term
    if LIN == 0
        cube = -(z).*(z.*z + zz.*zz);
%      % Linearized force (NOT YET IMPLEMENTED)
%      elseif LIN == 1
%          xxx
    end
    fo = a.*z + zz.*omega + cube;
    %%% HOPF (Langevin equation)
    z = z + dt*(fo + sum) + no;
end
z0 = z;
end
