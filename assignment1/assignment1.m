clear all 
close all 
clc 

%% Shubh Raval ME6103 Assignment 1 

%% Modified Rosenbrock Function
% f(x1,x2) = (2.8-x1)^2 + 100(x2-x1^2)^2

%% 2A: 
%For this problem I am chosing steepest descent which is a 
%method that minimizes a function by iteratively moving in the direction of the negative gradient. 
%It adjusts parameters step by step to reduce the function's value and find a local minimum.
%because it efficiently handles smooth, 
%convex-like regions and leverages gradient information for fast convergence. 

%% 2B Solve the problem: 

fun = @(x) (2.8 - x(1))^2 + 100 * (x(2) - x(1)^2)^2;

x0 = [-1.5, 2]; 

maxFuncEvals = 600;

error_delta = inf;


prev_fval = inf;

while error_delta > 0.001

    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
        'HessUpdate', 'steepdesc', ...
        'MaxFunctionEvaluations', maxFuncEvals, ... 
        'Display', 'iter'); 

    [x, fval, eflag] = fminunc(fun, x0, options);


    error_delta = abs(prev_fval - fval);

    if fval < prev_fval
        prev_fval = fval;
        x0 = x; 
    end

    maxFuncEvals = maxFuncEvals * 2;

    fprintf('New iteration: Error delta = %.6f, New x0 = [%.6f, %.6f], Max Func Evals = %d\n', ...
        error_delta, x0(1), x0(2), maxFuncEvals);
end

% Final results
disp('Validated Optimized Solution:');
disp(x0);
disp('Final Function Value:');
disp(prev_fval);

%% 2C Validation: 
function [f, grad, hessian] = modified_rosenbrock(x)
    f = (2.8 - x(1))^2 + 100 * (x(2) - x(1)^2)^2;

    if nargout > 1  % Gradient
        df_dx1 = -2 * (2.8 - x(1)) - 400 * x(1) * (x(2) - x(1)^2);
        df_dx2 = 200 * (x(2) - x(1)^2);
        grad = [df_dx1; df_dx2];
    end

    if nargout > 2  % Hessian
        d2f_dx1dx1 = 2 - 400 * (x(2) - 3 * x(1)^2);
        d2f_dx1dx2 = -400 * x(1);
        d2f_dx2dx1 = d2f_dx1dx2;
        d2f_dx2dx2 = 200;
        hessian = [d2f_dx1dx1, d2f_dx1dx2; d2f_dx2dx1, d2f_dx2dx2];
    end
end
disp('Performing Hessian Validation...');

options = optimoptions('fminunc', 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian', 'on', 'Display', 'iter');

[x_opt, fval, exitflag, output, grad, hessian] = fminunc(@modified_rosenbrock, x0, options);

% Display Results
disp('Optimized solution after Hessian validation:');
disp(x_opt);
disp('Function value at minimum:');
disp(fval);
disp('Hessian at minimum:');
disp(hessian);

% Check eigenvalues
eig_vals = eig(hessian);

if all(eig_vals > 0)
    disp('Validation Passed: Hessian is positive definite, solution is a local minimum.');
else
    disp('Hessian is not positive definite');
end


%% Modified Spring Function 
% f(x) = sin(x) + 0.05x^2

%% 2A
% For this problem I will once again be using the steepest descent method
% to reach convergence faster with a similar algorithim pushing the
% iterations to a final minima, this is a reasonable method since it is
% good for smooth functions where a descent over time should be feasible.

%% 2B Solve the Problem


fun = @(x) sin(x) + (0.05)*x^2 ; 

x0 = -1.0; 

maxFuncEvals = 600;


error_delta = inf;


prev_fval = inf;

while error_delta > 0.001

    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
        'HessUpdate', 'steepdesc', ... 
        'MaxFunctionEvaluations', maxFuncEvals, ... 
        'Display', 'iter');

    [x, fval, eflag, output] = fminunc(fun, x0, options);

    error_delta = abs(prev_fval - fval);


    if fval < prev_fval
        prev_fval = fval; 
        x0 = x; 
    end

    maxFuncEvals = maxFuncEvals * 2;

    fprintf('New iteration: Error delta = %.6f, New x0 = [%.6f], Max Func Evals = %d\n', ...
        error_delta, x0(1), maxFuncEvals);
end

% Final results
disp('Validated Optimized Solution:');
disp(x0);
disp('Final Function Value:');
disp(prev_fval);


%% 3 Unique Optimization Problem: 
% Cantilever beam deflection optimization problem: 

% Constants
E = 65e9; % Young's modulus (Pa)
F = 1100;   % Load at the tip (N)
target_deflection = 0.4; % Desired deflection (m)
rho = 2700; % Density of material (kg/m^3)
f_min = 4; % Minimum natural frequency (Hz)

x0 = [0.3, 0.05, 3];

lb = [0.2, 0.0254, 2];
ub = [0.5, 0.055, 5];

options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
x_opt = fmincon(@objective, x0, [], [], [], [], lb, ub, @(x) constraints(x, E, F, target_deflection, rho, f_min), options);
    
L_values = linspace(2, 5, 50);
deflections = arrayfun(@(L) (F * L^3) / (3 * E * ((x_opt(1) * x_opt(2)^3) / 12)), L_values);
frequencies = arrayfun(@(L) (1.875^2) * sqrt((E * ((x_opt(1) * x_opt(2)^3) / 12)) / (rho * (x_opt(1) * x_opt(2)) * L^4)), L_values);
volumes = arrayfun(@(L) x_opt(1) * x_opt(2) * L, L_values);

% Final Results
fprintf('Optimal Dimensions: Width = %.3f m, Thickness = %.3f m, Length = %.3f m\n', x_opt(1), x_opt(2), x_opt(3));

%% Validation
figure;
subplot(3,1,1);
plot(L_values, deflections, 'b-', 'LineWidth', 2);
hold on;
yline(target_deflection, 'r--');
xlabel('Beam Length (m)');
ylabel('Deflection (m)');
title('Deflection vs Length');

subplot(3,1,2);
plot(L_values, frequencies, 'g-', 'LineWidth', 2);
hold on;
yline(f_min, 'r--');
xlabel('Beam Length (m)');
ylabel('Natural Frequency (Hz)');
title('Frequency vs Length');

subplot(3,1,3);
plot(L_values, volumes, 'm-', 'LineWidth', 2);
xlabel('Beam Length (m)');
ylabel('Volume (m^3)');
title('Volume vs Length');

%% Seperated Objective and Constraint functions for clarity
function V = objective(x)
    V = x(1) * x(2) * x(3);
end

function [c, ceq] = constraints(x, E, F, target_deflection, rho, f_min)
    I = (x(1) * x(2)^3) / 12;
    deflection = (F * x(3)^3) / (3 * E * I);
    % ^^ Standard deflection eq
    A = x(1) * x(2);
    f = (1.875^2) * sqrt((E * I) / (rho * A * x(3)^4)); 
    % ^^ First natural frequency equation for a rectangular cross section of a CL Beam
    ceq = deflection - target_deflection;
    c = f_min - f;
end