function [ts, force, noise] = hopf_sim_1start(ahopf,omega,dsig,dt,Tmax,TR,NPARCELS,wC,sumC,LIN,z0)
% HOPF_SIM_1START Runs a single recorded Hopf simulation from initial state z0.
% Author: Juan Manuel Monti
%
% Integrates the Stuart-Landau network with Euler-Maruyama, sampling at every TR.
% Records x-component (BOLD proxy), deterministic force (negated per eq. 2.7 of
% Cugliandolo 1994), and noise divided by sqrt(dt) (continuous-time convention).
%
% Inputs:
%   ahopf      bifurcation parameter (fixed: -0.02)
%   omega      [NPARCELSx2] angular frequencies (rad/s); column 1 negated
%   dsig       noise amplitude pre-scaled by sqrt(dt)
%   dt         integration timestep (s)
%   Tmax       number of TR steps to record
%   TR         fMRI repetition time (s)
%   NPARCELS  number of brain parcels
%   wC         [NPARCELSxNPARCELS] effective connectivity (G*Cij)
%   sumC       [NPARCELSx2] row sums of wC replicated for x and y components
%   LIN        0 = cubic Hopf; 1 = linear (not yet implemented)
%   z0         [NPARCELSx2] initial state from hopf_sim_0init
%
% Outputs:
%   ts     [NPARCELSxTmax] x-component of z at TR intervals
%   force  [NPARCELSxTmax] negated deterministic force (Cugliandolo 1994 eq. 2.7)
%   noise  [NPARCELSxTmax] noise term divided by sqrt(dt)

%%% Start simulation from initialized z0
z = z0;

a = ahopf.*ones(NPARCELS,2);
xs = zeros(Tmax,NPARCELS);
ys = zeros(Tmax,NPARCELS);
force = zeros(Tmax,NPARCELS);
noise = zeros(Tmax,NPARCELS);

%%% Actual modeling (x = BOLD signal (Interpretation), y some other oscillation)
Tsim = (Tmax-1)*TR;
nn = 0;
for t = 0:dt:(Tsim)
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
        cube = -z.*(z.*z + zz.*zz);
%      % Linearized force (NOT YET IMPLEMENTED)
%      elseif LIN == 1
%          xxx
    end
    fo = a.*z + zz.*omega + cube;
    %%% HOPF (Langevin equation)
    z = z + dt*(fo + sum) + no;

    %%% Saves simulated TS, "force" and noise for this simulation
    if abs(mod(t,TR))<0.01
        nn = nn+1;
        xs(nn,:) = z(:,1)';
        ys(nn,:) = z(:,2)';
        % In order to keep the same definition as in Cugliandolo 1994
        % dv/dt = force[v](t) + sigma*eta(t)
        force(nn,:) = (fo(:,1)' + sum(:,1)');
        noise(nn,:) = no(:,1)';
    end
end
ts = xs';
force = -force'; % To be consistent with equation (2.7) cugliandolo 1994
% In the Euler–Maruyama method the noise is multiplied by sqrt(dt)
noise = noise'./sqrt(dt);
end
