clear 
close all 
clc

load = 50;
platform_diam = 0.09;
base_diam = 0.096; 
H = 0.03;



[mass, isValid, penalty, W_opt, T_opt, a_opt, b_opt, deflection_mm] = Deflection_Optim(load, H, platform_diam, base_diam);

disp(W_opt)
disp(T_opt)
disp(a_opt)
disp(b_opt)
disp(mass)