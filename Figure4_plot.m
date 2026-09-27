clear all
clc

% Figure 4: comparison between the full two-state model and the
% effective single-state model for R = 4, 5.66, and 8 um.

%% Load data
load('PDEresults.mat')
load('PDEresults_SingleState.mat')

% Data indices: R = 4, 5.66, and 8 um
id = [11 12 13];
Rlabel = {'4','5.7','8'};

%% Plot settings (matched to Figure3_plot.m)
lwBlack = 2;
lwRed   = 3;   % single-state lines slightly thicker
lwAxis  = 1;
fs      = 23;

Position = [100 100 900 700];
AxisLim  = [0 150 0 1];
XTick    = 0:50:150;
YTick    = 0:0.2:1;

% Base colors
blackBase = [0 0 0];
% Brighter red, closer to pure red while remaining slightly softened
redBase   = [235 45 35]/255;

% Apparent transparency levels, from light to dark as R increases.
% MATLAB line objects do not consistently support alpha transparency,
% so the colors are blended with white to give the same visual effect
% on a white background.
alphaLevel = [0.40 0.68 1.00];
white      = [1 1 1];

blackColor = zeros(3,3);
redColor   = zeros(3,3);
for i = 1:3
    blackColor(i,:) = alphaLevel(i)*blackBase + (1-alphaLevel(i))*white;
    redColor(i,:)   = alphaLevel(i)*redBase   + (1-alphaLevel(i))*white;
end

fig4 = figure(4);
clf(fig4)
set(fig4,'Position',Position,'Color','w')
hold on

h = gobjects(6,1);

for i = 1:3
    % Two-state model
    t_two  = PDEresults(id(i)).t;
    Ic_two = PDEresults(id(i)).Ic;

    % Single-state model
    t_single  = PDEresults_single(id(i)).t;
    Ic_single = PDEresults_single(id(i)).Ic;

    % Plot lighter to darker for increasing condensate radius
    h(2*i-1) = plot(t_two, Ic_two, '-', ...
        'LineWidth',lwBlack,'Color',blackColor(i,:), ...
        'DisplayName',[' two-state, {\it R} = ' Rlabel{i} ' \mum']);

    h(2*i) = plot(t_single, Ic_single, '--', ...
        'LineWidth',lwRed,'Color',redColor(i,:), ...
        'DisplayName',[' single-state, {\it R} = ' Rlabel{i} ' \mum']);
end

axis(AxisLim)
xticks(XTick)
yticks(YTick)

xlabel('{\it t} (s)');
ylabel('{\it I}_c({\itt})');

set(gca,'FontSize',fs,'LineWidth',lwAxis,'Box','on')

legend(h, 'Location','southeast', 'Box','off', ...
       'FontSize',fs-3, 'NumColumns',1);

% Optional export
% exportgraphics(fig4,'Figure4.pdf','ContentType','vector');
