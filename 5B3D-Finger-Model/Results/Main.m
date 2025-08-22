%% 5B3D-Finger-Model-Results: Main
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% This script helps analyze the results from the 5B3D-Model: Main

% clear all; close all; clc;
load('Results.mat');

%% Generate a video (.mp4) of the 
% rate = 10;
% name = 'Finger-Stable-Configurations.mp4';
% genVideo(param, N, t_1, t_2, t_3, rate, name);

%% Output Grasp Force Data for specific circular object diameter
% R = 101.6/1000/2; % input object diameter in m (4 in diameter)
% t_3 = linspace(deg2rad(-90), deg2rad(40), (90+40)*10+1)';
% [N_s, Nx, Ny, F] = circleCalc(param, R, t_3, 1);

%% Output Grasp Force Data for specific angles
% t_1 and t_2 are scalars
%% From Test Object in CAD
% t_3 = linspace(-0.8283, -0.6283, 201);
% t_1 = 1.0341;
% t_2 = 1.0363;
%% Another Test Object in CAD
% t_3 = linspace(0.2312, 0.2486, 2);
% t_1 = 1.5664;
% t_2 = 1.5446;
% [N_s, Nx, Ny, F] = ellipseCalc(param, t_1, t_2, t_3, 1);


%% Output the x and y-direction forces given specific enveloping force
% ref = [10 20]; % reference forces in N
% [Fx, Fy] = cartForces(F, Nx, Ny, ref)
