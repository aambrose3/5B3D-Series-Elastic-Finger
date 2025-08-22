%% Grasp_Force_Testing: GraspForceTestingAnalysis
%% Author: Alexander B. Ambrose
%% Date: 8-12-2025

%% Description:
%  This script is used to plot the differences between the modeled contact
%  force and the measureed contact force.
%  
%  Outputs:
%  N/A

clear all; close all; clc;

addpath '..\5B3D-Finger-Model'
addpath '..\5B3D-Finger-Model\Results'
load('Data/Matched_Data.mat'); 

start = 1;
for kk = 6:length(fin)
    if start == 1
        t1 = fin{kk}.t1;
        t2 = fin{kk}.t2;
        t3 = fin{kk}.t3;
        F = fin{kk}.F;
        Fx = rok{kk}.Fx;
        Fy = rok{kk}.Fz;
        start = 0;
    else
        t1 = [t1; fin{kk}.t1];
        t2 = [t2; fin{kk}.t2];
        t3 = [t3; fin{kk}.t3];
        F = [F; fin{kk}.F];
        Fx = [Fx; rok{kk}.Fx];
        Fy = [Fy; rok{kk}.Fz];
    end
end

%% Downsample to speed up computation
f = 1/(fin{1}.time(2) - fin{1}.time(1));
N = round(f/2); % 2 samples per second
t3_ = downsample(t3, N);
F_ = downsample(F, N);
F_x = downsample(Fx, N);
F_y = downsample(Fy, N);

%% Output Grasp Force Data for specific angles
% t_1 and t_2 need to be scalars and add uncertainties 
t_1 = [mean(t1)-2*std(t1); mean(t1)+2*std(t1)];
t_2 = [mean(t2)-2*std(t2); mean(t2)+2*std(t2)];
t_3 = linspace(min(t3)-deg2rad(0.125), max(t3)+deg2rad(0.125), 101)';
param = getParam();
[N_s, Nx, Ny, F_sim] = ellipseCalc(param, t_1, t_2, t_3, 1);
hold on
scatter(t3_, F_x, 50, '+b')
scatter(t3_, F_y, 50, '*r')
F_tot = sqrt(F_x.^2+F_y.^2);
scatter(t3_, F_tot, 50, 'ok')
ylabel('Grasping Contact Forces (N)', 'FontSize', 16)
xlabel('Actuator Angle (rad)', 'FontSize', 16)
title('Grasping Contact Forces Produced by One Finger', 'FontSize', 16)
legend('Predicted X', 'Predicted Y', 'Predicted Total',...
    'Measured X', 'Measured Y', 'Measured Total', ...
    'location', 'northwest', 'FontSize', 16)
axis([min(t3_) max(t3_) 0 40])


F_tot = sqrt(Fx.^2+Fy.^2);
x = 100*(F_tot-F)./F_tot;
figure
subplot(2, 1, 2)
boxplot(x, 'labels', 'N=5')
title('Total Enveloping Grasping Force Residuals', 'FontSize', 16)
ylabel('Percentage Error in Enveloping Grasp Force (%)', 'FontSize', 16)

subplot(2, 1, 1)
histogram(x)
title('Total Enveloping Grasping Force Error', 'FontSize', 16)
xlabel('Percentage Error in Enveloping Grasp Force (%)', 'FontSize', 16)
ylabel('Prevalence', 'FontSize', 16)

