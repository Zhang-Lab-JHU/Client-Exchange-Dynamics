clear all
clc

global Cbden Cbdil Cuden Cudil 
global Dbden Dbdil Duden Dudil  
global R0

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

koff = 0.3;     % unbinding rate k_off, s^-1
Taucon = 1/koff;

Au = Cuden/(Cuden + Cbden);
Ab = 1 - Au;

load('PDEresults.mat')

NR = length(PDEresults);
Tau1 = zeros(NR,1);
Tau2 = zeros(NR,1);
Tau3 = zeros(NR,1);
Radius = zeros(NR,1);

for n = 1:NR

    R0 = PDEresults(n).R0;
    t  = PDEresults(n).t;
    r  = PDEresults(n).r;
    cb = PDEresults(n).cb;
    cu = PDEresults(n).cu;
    Ic = PDEresults(n).Ic;

    [Taub,Tauu,Taue] = TauPrediction;

    ft = fittype('log(Au*exp(-t/Tau1)+(1-Au)*exp(-t/Tau2))', ...
                 'independent','t', ...
                 'problem','Au', ...
                 'coefficients',{'Tau1','Tau2'});
    fitObj = fit(t',log(1-Ic),ft,'problem',Au,'StartPoint',[Tauu,Taub]);
    Tau1(n) = fitObj.Tau1;
    Tau2(n) = fitObj.Tau2;

    ft = fittype('-t/Tau3', ...
                 'independent','t', ...
                 'coefficients','Tau3');
    fitObj = fit(t',log(1-Ic),ft,'StartPoint',Taue);
    Tau3(n) = fitObj.Tau3;

    Radius(n) = R0;
end

R0 = 1;
RR = 0.1:0.01:10;
[tau_b0, tau_u0, tau_e0] = TauPrediction;
tau_b   = tau_b0 * RR.^2;        % s
tau_u   = tau_u0 * RR.^2;        % s
tau_e   = tau_e0 * RR.^2;        % s
tau_c = (1/koff) * ones(size(RR)); % s
tau_slow = tau_b.*(tau_e+tau_c)./(tau_b+tau_c);
tau_fast = tau_u.*(tau_e+tau_c)./(tau_u+tau_c);


%%
id1 = 1:10;
id3 = 11:13;

lw1=1.5;
lw2=1.5;
lw3=1;
fs=23;
%Position=[100 100 800 800];
Position=[100 100 882.5 700];
Axis=[0.1 10 0.002 200];
XTick=[0.125 0.25 0.5 1 2 4 8];
YTick=[0.01 0.1 1 10 100];

fig2=figure(2);
set(fig2,'Position',Position)

color1=[0,114,189]/255;
color2=[217,83,25]/255;
color3=[78,167,46]/255;%[146,208,80]/255;
color4=[0,0,0]/255;

loglog(RR, tau_b,'--','LineWidth',lw2,'color','k'); hold on;
loglog(RR, tau_u,'-.','LineWidth',lw2,'color','k');
loglog(RR, tau_c,':','LineWidth',lw2,'color','k'); 

x_dot = logspace(log10(min(RR)),log10(max(RR)),100);
y_dot = interp1(RR,tau_e,x_dot,'pchip');
plot(x_dot,y_dot,'k.','MarkerSize',5);

loglog(RR, tau_b.*(tau_e+tau_c)./(tau_b+tau_c),'-','LineWidth',lw1,'color',color2);
loglog(RR, tau_u.*(tau_e+tau_c)./(tau_u+tau_c),'-','LineWidth',lw1,'color',color1);
scatter(Radius,[Tau2(id1);Tau3(id3)],400,'>',"filled",'MarkerFaceColor',color2,'MarkerEdgeColor',color2,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); hold on
scatter(Radius,[Tau1(id1);Tau3(id3)],400,'<',"filled",'MarkerFaceColor',color1,'MarkerEdgeColor',color1,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); hold on

axis(Axis)
xticks(XTick)
yticks(YTick)
% Labels and legend
xlabel('Condensate radius, {\it R} (\mum)');
ylabel('Timescale, {\tau} (s)');
% legend('\tau_b','\tau_u','\tau_{c}','\tau_{e}',...
%        '\tau_{slow}','\tau_{fast}','\tau_{slow}','\tau_{fast}',...
%        'Location','southeast','box','off');
set(gca,'FontSize',fs,'LineWidth',lw3)

function [Taub,Tauu,Taue] = TauPrediction
global Cbden Cbdil Cuden Cudil Dbden Dbdil Duden Dudil R0
Cden = Cbden+Cuden;
Cdil = Cbdil+Cudil;
Dcden = (Dbden*Cbden+Duden*Cuden)/Cden;
Dcdil = (Dbdil*Cbdil+Dudil*Cudil)/Cdil;
Taub = R0^2/(pi^2*Dbden)+Cbden*R0^2/(3*Cbdil*Dbdil);
Tauu = R0^2/(pi^2*Duden)+Cuden*R0^2/(3*Cudil*Dudil);
Taue = R0^2/(pi^2*Dcden)+Cden*R0^2/(3*Cdil*Dcdil);
end