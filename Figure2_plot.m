clear all
clc

global Cbden Cbdil Cuden Cudil 
global Dbden Dbdil Duden Dudil  
global R0 Rl

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

load('PDEresults.mat')

idxR = 5;   % Radius = 2.^(-3:0.5:3); entry 5 is R = 0.5 um
R0 = PDEresults(idxR).R0;
t  = PDEresults(idxR).t(:);
r  = PDEresults(idxR).r(:);
cb = PDEresults(idxR).cb;
cu = PDEresults(idxR).cu;
Ic = PDEresults(idxR).Ic(:);

Rl = 0.005;
[Cb_eq, Cu_eq, Db, Du] = Profiles(r);
[Taub, Tauu, Taue] = TauPrediction;
[tau_fast, tau_slow, Au] = ExtractTimescales(t, Ic, Tauu, Taub);
Ab = 1 - Au;
% Au = mean(Cu_eq(r<R0+Rl))/(mean(Cb_eq(r<R0+Rl))+mean(Cu_eq(r<R0+Rl)))

fprintf('R = %.3g um\n', R0);
fprintf('Fitted tau_fast = %.4g s, tau_slow = %.4g s\n', tau_fast, tau_slow);

% Time points for profiles
times_b = [0, tau_slow/4, tau_slow/2, tau_slow, 2*tau_slow];
times_u = [0, tau_fast/4, tau_fast/2, tau_fast, 2*tau_fast];
idx_b = TimeIndex(t, times_b);
idx_u = TimeIndex(t, times_u);

%% Plot settings
% For panels a and b: unbound/pore-space clients are blue; bound/scaffold-associated clients are orange.
color_b = [217, 83, 25]/255;
color_u = [0, 114, 189]/255;
color_fit = [217, 83, 25]/255;
profile_colors = [0.00 0.00 0.00;
                  0.05 0.45 0.15;
                  0.20 0.65 0.20;
                  0.35 0.75 0.35;
                  0.10 0.70 0.10];
profile_styles = {'-', '--', ':', '-.', '--'};

lw = 1.4;
lwFit = 1.6;
lwAxis = 1.0;
fs = 16;
fsLabel = 17;
fsPanel = 18;
fsLegend = 13.5;

cell = [100 0 650 800];

% Pixel layout requested by user.
% Each of a--d occupies a 250 x 250 cell.  Horizontal spacing is 50 px.
% Vertical spacing between rows a/b and c/d is 25 px, and between c/d and e is 25 px.
% Panel e occupies a 550 x 250 cell.

gap = 30;
gapAB = 50;
gapAC = 10;
sizeX = (cell(3)-gapAB)/2;
sizeY = (cell(4)-2*gapAC)/3;
h = 30;
d = 10;
cellA = [gap   gap+2*(sizeY+gapAC)+h-d sizeX sizeY+h];
cellB = [gap+sizeX+gapAB gap+2*(sizeY+gapAC)+h-d sizeX sizeY+h];
cellC = [gap   gap+sizeY+gapAC-d sizeX sizeY+h];
cellD = [gap+sizeX+gapAB gap+sizeY+gapAC-d sizeX sizeY+h];
cellE = [gap     gap sizeX*2+gapAB sizeY-d];

cell(3) = cell(3)+gap;
cell(4) = cell(4)+gap+2*h;

fig = figure(2);
clf(fig)
set(fig, 'Units', 'pixels', 'Position', cell, 'Color', 'w')

% Axes inside each cell. These margins leave room for y-labels, x-labels,
% and panel letters outside the box while preserving the requested cell sizes.
posA = cellAxes(cellA, 54, 48, 15, 28);
posB = cellAxes(cellB, 54, 48, 15, 28);
posC = cellAxes(cellC, 54, 48, 15, 28);
posD = cellAxes(cellD, 54, 48, 15, 28);
posE = cellAxes(cellE, 54, 48, 15, 28);

% Broken-x setup. Left segment spans 0--1 um; right segment spans 24.5--25 um.
% The right segment is displayed with half the width of the left segment,
% matching the physical interval lengths (1 um vs 0.5 um).
xLeft  = [0, 1.2];
xRight = [24.7, 25.0];
mapOpt.leftEnd = xLeft(2);
mapOpt.rightStart = xRight(1);
mapOpt.rightEnd = xRight(2);
mapOpt.gap = 0.090;  % reduced spacing between the two break marks
mapOpt.rightWidth = xRight(2) - xRight(1);  % equal scale with left segment
mapOpt.xMax = mapOpt.leftEnd + mapOpt.gap + mapOpt.rightWidth;

xTicksPhysical = [0, 0.5, 1, 25.0];
xTicksMapped = mapBrokenX(xTicksPhysical, mapOpt);
xTickLabels = {'0', '0.5', '1', '25'};

% a. Equilibrium concentration profiles
axA = axes(fig, 'Units', 'pixels', 'Position', posA); hold(axA, 'on')
plotBroken(axA, r, Cb_eq, mapOpt, '-', color_b, lwFit)
plotBroken(axA, r, Cu_eq, mapOpt, '-', color_u, lwFit)
formatBrokenPanel(axA, mapOpt, [0.05 200], 'log', xTicksMapped, xTickLabels, fs, lwAxis)
ylabel(axA, 'Equilibrium conc. (\muM)', 'FontSize', fsLabel)
legend(axA, {'{\it c_b}^{eq}', '{\it c_u}^{eq}'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', fsLegend)
xlabel(axA, '{\it r} (\mum)', 'FontSize', fsLabel)
%panelLabel(axA, 'a', fsPanel)

% b. Diffusion coefficient profiles
axB = axes(fig, 'Units', 'pixels', 'Position', posB); hold(axB, 'on')
plotBroken(axB, r, Db, mapOpt, '-', color_b, lwFit)
plotBroken(axB, r, Du, mapOpt, '-', color_u, lwFit)
formatBrokenPanel(axB, mapOpt, [0.05 200], 'log', xTicksMapped, xTickLabels, fs, lwAxis)
ylabel(axB, 'Diffusion coef. (\mum^2/s)', 'FontSize', fsLabel)
legend(axB, {'{\it D_b}', '{\it D_u}'}, 'Location', 'southeast', 'Box', 'off', 'FontSize', fsLegend)
xlabel(axB, '{\it r} (\mum)', 'FontSize', fsLabel)
%panelLabel(axB, 'b', fsPanel)

% c. Bound-client concentration profiles
axC = axes(fig, 'Units', 'pixels', 'Position', posC); hold(axC, 'on')
for k = 1:numel(idx_b)
    plotBroken(axC, r, cb(idx_b(k),:), mapOpt, profile_styles{k}, profile_colors(k,:), lw)
end
cbMax = 1.05*max(cb(idx_b(1),:));
formatBrokenPanel(axC, mapOpt, [0 cbMax], 'linear', xTicksMapped, xTickLabels, fs, lwAxis)
ylabel(axC, '{\it c_b}({\itr},{\it t}) (\muM)', 'FontSize', fsLabel)
legend(axC, {'{\it t} = 0', '{\it t} = \tau_{slow}/4', '{\it t} = \tau_{slow}/2', ...
             '{\it t} = \tau_{slow}', '{\it t} = 2\tau_{slow}'}, ...
             'Location', 'northeast', 'Box', 'off', 'FontSize', fsLegend)
xlabel(axC, '{\it r} (\mum)', 'FontSize', fsLabel)
%panelLabel(axC, 'c', fsPanel)
yticks(0:25:100)

% d. Unbound-client concentration profiles
axD = axes(fig, 'Units', 'pixels', 'Position', posD); hold(axD, 'on')
for k = 1:numel(idx_u)
    plotBroken(axD, r, cu(idx_u(k),:), mapOpt, profile_styles{k}, profile_colors(k,:), lw)
end
cuMax = 1.05*max(cu(idx_u(1),:));
formatBrokenPanel(axD, mapOpt, [0 cuMax], 'linear', xTicksMapped, xTickLabels, fs, lwAxis)
ylabel(axD, '{\it c_u}({\itr},{\it t}) (\muM)', 'FontSize', fsLabel)
legend(axD, {'{\it t} = 0', '{\it t} = \tau_{fast}/4', '{\it t} = \tau_{fast}/2', ...
             '{\it t} = \tau_{fast}', '{\it t} = 2\tau_{fast}'}, ...
             'Location', 'northeast', 'Box', 'off', 'FontSize', fsLegend)
xlabel(axD, '{\it r} (\mum)', 'FontSize', fsLabel)
%panelLabel(axD, 'd', fsPanel)
yticks(0:25:100)

% Add the broken-axis graphics after the axes have been created so the white
% gap is drawn above the axis box.
addBreakGraphics(fig, [axA axB axC axD], mapOpt, 1.05)

% e. Full-FRAP recovery curve with short-time inset
axE = axes(fig, 'Units', 'pixels', 'Position', posE);
hold(axE, 'on')
plot(axE, t, Ic, 'k-', 'LineWidth', lwFit)
I_fit = 1 - (Au*exp(-t/tau_fast) + Ab*exp(-t/tau_slow));
plot(axE, t, I_fit, '--', 'Color', color_fit, 'LineWidth', lwFit*1.5)
plot(axE, tau_slow, 0, '.');
box(axE, 'on')
set(axE, 'FontSize', fs, 'LineWidth', lwAxis, 'TickDir', 'in')
xlabel(axE, '{\it t} (s)', 'FontSize', fsLabel)
ylabel(axE, '{\it I_c}({\itt})', 'FontSize', fsLabel)
legend(axE, {'{ simulation}', '{ fit}'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', fsLegend)
%xlim(axE, [0, max(t)])
ylim(axE, [0, 1.02])
yticks(0:0.25:1)
%panelLabel(axE, 'e', fsPanel)

% Inset: short-time recovery on the fast timescale; place at lower right of panel e.
insetW = 200;
insetH = 100;
insetPos = [posE(1)+posE(3)-insetW-24, posE(2)+40, insetW, insetH];
axInset = axes(fig, 'Units', 'pixels', 'Position', insetPos);
hold(axInset, 'on')
plot(axInset, t, Ic, 'k-', 'LineWidth', lwFit)
plot(axInset, t, I_fit, '--', 'Color', color_fit, 'LineWidth', lwFit*1.5)
box(axInset, 'on')
set(axInset, 'FontSize', fsLegend, 'LineWidth', lwAxis, 'TickDir', 'in')
xlim(axInset, [0, 3*tau_fast])
ylim(axInset, [0, 1.05*max(Ic(t <= 3*tau_fast))])
yticks(0:0.25:1)
plot(axInset, tau_fast, 0, '.');
%xlabel(axInset, '{\it t} (s)', 'FontSize', fsLegend)
%ylabel(axInset, '{\it I_c}({\itt})', 'FontSize', fsLegend)





%%
function [Cb, Cu, Db, Du] = Profiles(r)
    global R0 Rl Cbden Cbdil Cuden Cudil Dbden Dbdil Duden Dudil
    z = (r - R0)/Rl;
    tanhz = tanh(z);
    Cb = 0.5*(Cbdil - Cbden).*(tanhz + 1) + Cbden;
    Cu = 0.5*(Cudil - Cuden).*(tanhz + 1) + Cuden;
    Db = 0.5*(Dbdil - Dbden).*(tanhz + 1) + Dbden;
    Du = 0.5*(Dudil - Duden).*(tanhz + 1) + Duden;
end

function [Taub, Tauu, Taue] = TauPrediction
    global Cbden Cbdil Cuden Cudil Dbden Dbdil Duden Dudil R0
    Cden = Cbden + Cuden;
    Cdil = Cbdil + Cudil;
    Dcden = (Dbden*Cbden + Duden*Cuden)/Cden;
    Dcdil = (Dbdil*Cbdil + Dudil*Cudil)/Cdil;
    Taub = R0^2/(pi^2*Dbden) + Cbden*R0^2/(3*Cbdil*Dbdil);
    Tauu = R0^2/(pi^2*Duden) + Cuden*R0^2/(3*Cudil*Dudil);
    Taue = R0^2/(pi^2*Dcden) + Cden*R0^2/(3*Cdil*Dcdil);
end

function [Tau_fast, Tau_slow, Au] = ExtractTimescales(t, Ic, Tauu, Taub)
    ft = fittype('log(Au*exp(-t/Tau1)+(1-Au)*exp(-t/Tau2))', ...
                'independent','t', ...
                'coefficients',{'Tau1','Tau2','Au'});
    fitObj = fit(t,log(1-Ic), ft, 'StartPoint', [Tauu, Taub, 0.5])
    Tau_fast = fitObj.Tau1;
    Tau_slow = fitObj.Tau2;
    Au = fitObj.Au;
end

function idx = TimeIndex(t, timesWanted)
    idx = zeros(size(timesWanted));
    for i = 1:numel(timesWanted)
        [~, idx(i)] = min(abs(t - timesWanted(i)));
    end
end

%% ---------------- Helper functions ----------------
function axPos = cellAxes(cellPos, leftPad, bottomPad, rightPad, topPad)
    axPos = [cellPos(1)+leftPad, cellPos(2)+bottomPad, ...
             cellPos(3)-leftPad-rightPad, cellPos(4)-bottomPad-topPad];
end

function xm = mapBrokenX(x, opt)
    xm = nan(size(x));
    leftMask = x <= opt.leftEnd;
    rightMask = x >= opt.rightStart;
    xm(leftMask) = x(leftMask);
    xm(rightMask) = opt.leftEnd + opt.gap + opt.rightWidth * ...
                    (x(rightMask) - opt.rightStart) / (opt.rightEnd - opt.rightStart);
end

function plotBroken(ax, x, y, opt, ls, color, lw)
    leftMask = x <= opt.leftEnd;
    rightMask = x >= opt.rightStart;

    % The curve is drawn in two pieces because of the broken x-axis.
    % Only the left piece should be visible to legend; otherwise MATLAB
    % counts the right piece as a duplicate legend entry.
    plot(ax, mapBrokenX(x(leftMask), opt), y(leftMask), ls, ...
         'Color', color, 'LineWidth', lw)
    plot(ax, mapBrokenX(x(rightMask), opt), y(rightMask), ls, ...
         'Color', color, 'LineWidth', lw, 'HandleVisibility', 'off')
end

function formatBrokenPanel(ax, opt, yLim, yScale, xTicksMapped, xTickLabels, fs, lwAxis)
    box(ax, 'on')
    set(ax, 'FontSize', fs, 'LineWidth', lwAxis, 'TickDir', 'in', ...
            'YScale', yScale, 'XLim', [0 opt.xMax], 'YLim', yLim, ...
            'XTick', xTicksMapped, 'XTickLabel', xTickLabels, ...
            'TickLength', [0.018 0.018], 'Layer', 'bottom')

    if strcmp(yScale, 'log')
        set(ax, 'YTick', [0.1 1 10 100])
    else
        set(ax, 'YTick', [0 50 100])
    end
end

function addBreakGraphics(fig, axList, opt, lwAxis)
    % Draw white gaps and black slash marks in figure-normalized coordinates.
    % This is more reliable than trying to delete part of the MATLAB box line.
    for ax = axList
        x1 = opt.leftEnd;
        x2 = opt.leftEnd + opt.gap;
        [xf1, yBot] = dataToFigNorm(fig, ax, x1, 0);
        [xf2, yTop] = dataToFigNorm(fig, ax, x2, 1);
        pos = ax.Position;
        figPos = fig.Position;
        posN = [pos(1)/figPos(3), pos(2)/figPos(4), pos(3)/figPos(3), pos(4)/figPos(4)];
        yB = posN(2);
        yT = posN(2) + posN(4);

        % White out the top and bottom box line only between the two slashes.
        annotation(fig, 'line', [xf1 xf2], [yB yB], 'Color', 'w', 'LineWidth', 3.2);
        annotation(fig, 'line', [xf1 xf2], [yT yT], 'Color', 'w', 'LineWidth', 3.2);

        % Compact diagonal marks, placed at the two edges of the blank gap.
        dx = 0.012;
        dy = 0.018;
        for xf = [xf1 xf2]
            annotation(fig, 'line', [xf-dx/2 xf+dx/2], [yB-dy/2 yB+dy/2], ...
                       'Color', 'k', 'LineWidth', lwAxis);
            annotation(fig, 'line', [xf-dx/2 xf+dx/2], [yT-dy/2 yT+dy/2], ...
                       'Color', 'k', 'LineWidth', lwAxis);
        end
    end
end

function [xf, yf] = dataToFigNorm(fig, ax, xData, yNorm)
    % Convert a data x coordinate and normalized y coordinate within the axes
    % to normalized figure coordinates. This is used only for annotations.
    oldUnitsAx = ax.Units;
    oldUnitsFig = fig.Units;
    ax.Units = 'pixels';
    fig.Units = 'pixels';
    pos = ax.Position;
    figPos = fig.Position;
    xl = ax.XLim;
    xf = (pos(1) + (xData - xl(1))/diff(xl)*pos(3)) / figPos(3);
    yf = (pos(2) + yNorm*pos(4)) / figPos(4);
    ax.Units = oldUnitsAx;
    fig.Units = oldUnitsFig;
end

function panelLabel(ax, txt, fsPanel)
    text(ax, -0.30, 1.10, txt, 'Units', 'normalized', 'FontSize', fsPanel, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'bottom', 'Clipping', 'off')
end
