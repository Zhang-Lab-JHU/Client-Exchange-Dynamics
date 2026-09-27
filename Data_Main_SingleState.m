function Data_Main_SingleState
% Effective single-state client exchange in a spherical condensate
% using MATLAB pdepe.
%
% This solves Eq. (22):
%   dc/dt = div{ Dc [grad(c) - c grad(c_eq)/c_eq] }
%
% with
%   c_eq(r) = Cb(r) + Cu(r)
%   Dc(r)   = [Db(r) Cb(r) + Du(r) Cu(r)] / c_eq(r)
%
% Parameters and profiles match Data_Reviewer.m for direct comparison.

clc

%% Parameters

% Equilibrium concentrations, uM
Cbden = 100;
Cbdil = 0.1;
Cuden = 100;
Cudil = 1;

% Diffusion coefficients, um^2/s
Dbden = 0.1;
Dbdil = 10;
Duden = 5;
Dudil = 100;

Rl = 0.005;      % interface width parameter
koff = 0.3;      % s^-1, used only to choose the same simulation time window
Taucon = 1/koff;

Radius = 2.^(-3:0.5:3);
NR = length(Radius);

PDEresults_single = struct();

for n = 1:NR

    n
    R0 = Radius(n);
    
    Rm = R0 - 2*Rl;
    Rp = R0 + 2*Rl;
    Rout = 50*R0;

    r = rMesh;

    [Taub,Tauu,Taue] = TauPrediction; %#ok<ASGLU>

    % Same time-window choice as in Data_Reviewer.m
    tEnd = (Taub*(Taue + Taucon)/(Taub + Taucon))*3;
    t = [0, logspace(-5,log10(tEnd),1000)];

    m = 2; % 3D spherical symmetry

    tic
    sol = pdepe(m, @pde, @pdeic, @pdebc, r, t);
    toc

    c = sol(:,:,1);
    Ic = Recovery(r,c);

    PDEresults_single(n).R0 = R0;
    PDEresults_single(n).t = t;
    PDEresults_single(n).r = r;
    PDEresults_single(n).c = c;
    PDEresults_single(n).Ic = Ic;

end

save('PDEresults_SingleState.mat','PDEresults_single','-v7.3');

%% ------------------------------------
% Make nonuniform radial mesh
function r = rMesh

    dr1 = Rl/2;
    dr2 = Rl/5;

    r1 = 0:dr1:Rm;
    r2 = Rm:dr2:Rp;
    r3 = Rp:dr1:Rout;

    r = unique([r1, r2, r3]);

end


%% ------------------------------------
% Effective single-state PDE, Eq. (22)
function [ccoeff,f,s] = pde(r,t,u,dudr) 

    z = (r - R0)/Rl;

    tanhz = tanh(z);
    sech2 = 1./cosh(z).^2;

    Cb = 0.5*(Cbdil - Cbden)*(tanhz + 1) + Cbden;
    Cu = 0.5*(Cudil - Cuden)*(tanhz + 1) + Cuden;

    Db = 0.5*(Dbdil - Dbden)*(tanhz + 1) + Dbden;
    Du = 0.5*(Dudil - Duden)*(tanhz + 1) + Duden;

    Ceq_local = Cb + Cu;

    dCbdr = 0.5*(Cbdil - Cbden)*sech2/Rl;
    dCudr = 0.5*(Cudil - Cuden)*sech2/Rl;
    dCeqdr = dCbdr + dCudr;

    % Eq. (23)
    Dc = (Db*Cb + Du*Cu)/Ceq_local;

    ccoeff = 1;
    f = Dc*(dudr - u*dCeqdr/Ceq_local);
    s = 0;

end


%% ------------------------------------
% Initial condition
function u0 = pdeic(r)

    u0 = Ceq(r).*(r <= Rp);

end


%% ------------------------------------
% No-flux boundary conditions
function [pl,ql,pr,qr] = pdebc(rl,ul,rr,ur,t) %#ok<INUSD>

    pl = 0;
    ql = 1;

    pr = 0;
    qr = 1;

end


%% ------------------------------------
% Equilibrium profiles
function y = Cb(r)

    y = 0.5*(Cbdil - Cbden).*(tanh((r - R0)/Rl) + 1) + Cbden;

end

function y = Cu(r)

    y = 0.5*(Cudil - Cuden).*(tanh((r - R0)/Rl) + 1) + Cuden;

end

function y = Ceq(r)

    y = Cb(r) + Cu(r);

end


%% ------------------------------------
% Analytical timescales
function [Taub,Tauu,Taue] = TauPrediction

    Cden = Cbden + Cuden;
    Cdil = Cbdil + Cudil;

    Dcden = (Dbden*Cbden + Duden*Cuden)/Cden;
    Dcdil = (Dbdil*Cbdil + Dudil*Cudil)/Cdil;

    Taub = R0^2/(pi^2*Dbden) + Cbden*R0^2/(3*Cbdil*Dbdil);
    Tauu = R0^2/(pi^2*Duden) + Cuden*R0^2/(3*Cudil*Dudil);
    Taue = R0^2/(pi^2*Dcden) + Cden*R0^2/(3*Cdil*Dcdil);

end


%% ------------------------------------
% Compute full-FRAP recovery curve
function Ic = Recovery(r,c)

    idx = r <= Rp;
    rr = r(idx);

    denom = trapz(rr, Ceq(rr).*rr.^2);

    Nt = size(c,1);
    Ic = zeros(Nt,1);

    for it = 1:Nt
        numer = trapz(rr, c(it,idx).*rr.^2);
        Ic(it) = 1 - numer/denom;
    end

end

end
