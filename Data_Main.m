function Data_Main
% Client exchange in a spherical condensate using MATLAB pdepe
%
% Variables:
%   u(1) bleached bound-client concentration
%   u(2) bleached unbound-client concentration

clc

% Parameters

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

Rl = 0.005;     % interface width parameter, um = 5 nm
koff = 0.3;     % unbinding rate k_off, s^-1
Taucon = 1/koff;

Radius = 2.^(-3:0.5:3);
NR = length(Radius);

PDEresults = struct();

for n = 1:NR

    n
    R0 = Radius(n);      % condensate radius, um
    Rm = R0 - 2*Rl;      % interface lower bound, um
    Rp = R0 + 2*Rl;      % interface upper bound, um
    Rout = 50*R0;        % outer boundary, um

    r = rMesh;

    [Taub,Tauu,Taue] = TauPrediction;  % Analytical timescales, s
    tEnd = (Taub*(Taue + Taucon)/(Taub + Taucon))*3;   % s
    t = [0, logspace(-5,log10(tEnd),1000)];

    m = 2; % symmetry constant, 3D spherical coordinates
    tic
    sol = pdepe(m, @pde, @pdeic, @pdebc, r, t);
    toc

    cb = sol(:,:,1);
    cu = sol(:,:,2);
    Ic = Recovery(r, cb, cu);

    PDEresults(n).R0 = R0;
    PDEresults(n).t = t;
    PDEresults(n).r = r;
    PDEresults(n).cb = cb;
    PDEresults(n).cu = cu;
    PDEresults(n).Ic = Ic;

end

save('PDEresults.mat','PDEresults','-v7.3');

%------------------------------------
% Make nonuniform radial mesh
function r = rMesh

    dr1 = Rl/2;
    dr2 = Rl/5;

    r1 = 0:dr1:Rm;
    r2 = Rm:dr2:Rp;
    r3 = Rp:dr1:Rout;

    r = unique([r1, r2, r3]);

end

%------------------------------------
% Reaction-diffusion PDE
function [c,f,s] = pde(r,t,u,dudr)

    z = (r - R0)/Rl;

    tanhz = tanh(z);
    sech2 = 1./cosh(z).^2;

    Cb = 0.5*(Cbdil - Cbden)*(tanhz + 1) + Cbden;
    Cu = 0.5*(Cudil - Cuden)*(tanhz + 1) + Cuden;

    Db = 0.5*(Dbdil - Dbden)*(tanhz + 1) + Dbden;
    Du = 0.5*(Dudil - Duden)*(tanhz + 1) + Duden;

    dCbdr = 0.5*(Cbdil - Cbden)*sech2/Rl;
    dCudr = 0.5*(Cudil - Cuden)*sech2/Rl;

    c = [1; 1];

    f = [Db*(dudr(1) - u(1)*dCbdr/Cb); Du*(dudr(2) - u(2)*dCudr/Cu)];

    F = koff*(Cb*u(2)/Cu - u(1));

    s = [F; -F];

end

%------------------------------------
% Initial conditions
function u0 = pdeic(r)

    u0 = [Cb(r); Cu(r)].*(r <= Rp);

end

%------------------------------------
% Boundary conditions
function [pl,ql,pr,qr] = pdebc(rl,ul,rr,ur,t)

    pl = [0; 0]; % ignored by solver since m = 2
    ql = [1; 1]; % ignored by solver since m = 2

    pr = [0; 0];
    qr = [1; 1];

end

%------------------------------------
% Equilibrium profiles and derivatives
function y = Cb(r)

    y = 0.5*(Cbdil - Cbden).*(tanh((r - R0)/Rl) + 1) + Cbden;

end

function y = Cu(r)

    y = 0.5*(Cudil - Cuden).*(tanh((r - R0)/Rl) + 1) + Cuden;

end

%------------------------------------
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

%------------------------------------
% Compute full-FRAP recovery curve
function Ic = Recovery(r,cb,cu)

    idx = r <= Rp;
    rr = r(idx);

    denom = trapz(rr, (Cb(rr) + Cu(rr)).*rr.^2);

    Nt = size(cb,1);
    Ic = zeros(Nt,1);

    for it = 1:Nt
        c_bleached = cb(it,idx) + cu(it,idx);
        numer = trapz(rr, c_bleached.*rr.^2);
        Ic(it) = 1 - numer/denom;
    end

end

end