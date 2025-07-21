% Minimal_NSGA2.m
% Most basic version possible
clear; close all; clc;

% Parameters
params.load = 9.81;
params.platform_diam = 0.10;
params.base_rad = 0.16;
params.h = 0.20;
params.theta_deg = 15;
params.actuator_masses = [41, 90, 107, 120, 122, 135, 195, 207.8, 230, 248, 248, 260, 282, 300, 450, 450, 570, 4000];
params.debug = false;

% Initialize a single point
W = 0.01;
T = 0.01;
actuator_idx = 1;

% Evaluate it
x = [W, T, actuator_idx];
fprintf('Evaluating point W=%f, T=%f, actuator=%d\n', W, T, actuator_idx);
F = SPM_objectives(x, params);
fprintf('Result: Mass = %f, GCI = %f\n', F(1), -F(2));

% Done!
fprintf('Test completed successfully\n');