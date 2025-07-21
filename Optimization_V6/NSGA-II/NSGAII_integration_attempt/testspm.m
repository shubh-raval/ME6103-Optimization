% test_SPM.m
clear; clc;
% Define sample decision vector and parameter structure
x = [0.01, 0.01, 1]; % W_init, T_init, actuator selection index
params.load = 9.81;
params.h = 0.15;
params.theta_deg = 5;
params.platform_diam = 0.14;
params.base_rad = 0.12;
params.actuator_masses = [5, 7, 10];
params.debug = true; % set to true to see debug prints

% Call SPM_objectives
f = SPM_objectives(x, params);
disp(f);