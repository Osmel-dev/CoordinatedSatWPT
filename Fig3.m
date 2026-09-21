% This code produces Fig. 3 of the following paper

%O. M. Rosabal, A. Azarbahram, M. Ashraf, M. Shehab, A. B. Khattak, O. L. A. López, and M.-S. Alouini, “Power from space: Coordinated satellite charging
% for off-grid wireless systems,” arXiv preprint arXiv:2608.25589, 2026.

% Keep this and Simulation_func_PL.m and pointing_error_rv.m files in the
% same folder and then run this file.

% Written by: Mateen Ashraf

clear all;

close all;

clc;

r_aper = 10 * 10^(-2); % receiver aperture size

D = 10^(-2) * [1:1:500];

w_0 = D/2; % transmit beamwidth at the transmitter

lambda_light = 1550*10^(-9); % first laser wavelength

P_loss = zeros(1,length(w_0));

Prop_dist_all = [100, 1000, 10000, 100000]; % propagation distances (100 Km, 1000 Km,..)

PL_all = zeros(2*length(Prop_dist_all), length(D)); % variable for storing pointing loss

% storing the results for first laser wavelength

for p=1:length(Prop_dist_all)

    Prop_dist = Prop_dist_all(p);

    for k = 1:length(w_0)

        P_loss(k) = Simulation_func_PL(Prop_dist, r_aper, w_0(k), lambda_light); % calls the function to give pointing loss

    end

    PL_all(p,:) = -10*log10(P_loss); % storing the results for first laser wavelength

end

lambda_light = 660*10^(-9); % second laser wavelength

% storing the results for second laser wavelength

for p=1:length(Prop_dist_all)

    Prop_dist = Prop_dist_all(p);

    for k = 1:length(w_0)

        P_loss(k) = Simulation_func_PL(Prop_dist, r_aper, w_0(k), lambda_light); % Function for calculating pointing loss

    end

    PL_all(p+length(Prop_dist_all),:) = -10*log10(P_loss); % storing the results for second laser wavelength

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following is just for plotting the pointing loss

Q = 50; % used for reducing the number of markers later

figure

set(gcf, 'Units', 'centimeters');  

axesFontSize = 14;                

legendFontSize = 12;              

afFigurePosition = [2 7 19 12]; 

set(gcf, 'Position', afFigurePosition,'PaperSize',[19 11],'PaperPositionMode','auto');  

hold on;

p1 = plot(D, PL_all(1,:), '-','Color',[.635, .078, .184],'LineWidth',1.5);  

p2 = plot(D, PL_all(2,:),    '-', 'Color',[.466, .674, .188],'LineWidth',1.5);

p3 = plot(D, PL_all(4,:), '-', 'Color',[0,.447, .7410],'LineWidth',1.5);

p4 = plot(D, PL_all(5,:),    'r--','Color',[.635, .078, .184], 'LineWidth',1.5);

p5 = plot(D, PL_all(6,:) , '--', 'Color',[.466, .674, .188],'LineWidth',1.5);

p6 = plot(D, PL_all(8,:),    '--', 'Color',[0,.447, .7410], 'LineWidth',1.5);

p7 = plot(D, PL_all(3,:),    'k-', 'LineWidth',1.5);

p8 = plot(D, PL_all(7,:),    'k--', 'LineWidth',1.5);

ax = gca;

ax.TickLabelInterpreter = 'latex';

ax.FontSize = 14;

set(p1, 'Marker','+','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p2, 'Marker','o','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p3, 'Marker','^','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p4, 'Marker','+','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p5, 'Marker','o','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p6, 'Marker','^','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p7, 'Marker','s','MarkerIndices',1:Q:length(D), MarkerSize = 10);

set(p8, 'Marker','s','MarkerIndices',1:Q:length(D), MarkerSize = 10);

sel_res{1} = 100;

sel_res{2} = 1000;

sel_res{3} = 10000;

sel_res{4} = 100000;

legend_str = cell(1,3);

k = 1;

for i = 1:4

    legend_str{k}   = sprintf('$%d$ km', sel_res{i});

    k = k+1;

end

grid on

box on

hold on

h1 = plot(nan,nan,'k-','LineWidth',1.5);  

h2 = plot(nan,nan,'k--','LineWidth',1.5);

lgd1 = legend([h1 h2],{'1550 nm','660 nm'},'Interpreter','latex','FontSize',14);

lgd1.Title.String = 'laser wavelength';

lgd1.AutoUpdate = 'off';

hSq = plot(nan,nan,'+','Color',[.635, .078, .184],'LineStyle','none',MarkerSize = 10);

hUp = plot(nan,nan,'o','Color',[.466, .674, .188],'LineStyle','none',MarkerSize = 10);

hCr = plot(nan,nan,'^', 'Color',[0,.447, .7410],'LineStyle','none',MarkerSize = 8);

hLM = plot(nan,nan,'ks','LineStyle','none',MarkerSize = 10);

lgd1.Position = [0.485 0.77 .15 0.10];

xlabel('transmit lens diameter (m)','Interpreter','latex', 'FontSize',14)

ylabel('pointing loss (dB)','Interpreter','latex', 'FontSize',14)

ax1 = gca;

ax2 = axes('Position',ax1.Position,'Visible','off');

lgd2 = legend(ax2,[hSq, hUp, hLM, hCr],legend_str,'Interpreter','latex','FontSize',14);

lgd2.Title.String = 'Link distance';

lgd2.Position = [0.665 0.67 0.2 0.20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%