% Run through a designed sequence for the SPM model
% It helps runtime to manually open the SPM Simulink model first.
% clear; close all; clc
addpath("Actuator","Kinematics","NSGA-II","SPM","Deflection_optim");
bdclose('SPM_Model');    % Close the model
load_system('SPM_Model'); % Re-open it clean
set_param('SPM_Model', 'FastRestart', 'off'); % Disable fast restart

%% Create the design sequence
N = 400;
p = haltonset(3);
z_n = net(p, N+1);
x_n = 1 - z_n(2:end,:);

% x := [h (m), theta (deg), diam_base (m)]
lb = [0.01, 30, .095];
ub = [.1, 150, .125];
x = x_n.*(ub - lb) + lb;

%% Iterate through calculating the values

% Define constants
diam_plat = .09; %m

T = nan(N,7);

for i = 1:N
    fprintf('i = %i -- Inputs: h = %.4f m, theta = %.2f deg, diam_base = %.4f.\n',i, x(i,1),x(i,2), x(i,3));

    diam_base = x(i,3); %m
    try
    % Run the SPM and deflection models to compute SPM mass
    [Torque, Force] = SPM_Simulation(x(i,1),x(i,2),diam_base,diam_plat);
    [m_act, ~] = actuator_selection(Torque);
    [m_link, isValid, penalty, ~,~,~,~, deflection_mm] = Deflection_Optim(max(Force.Data(end,:)), x(i,1), diam_plat, diam_base, false);
    
    % Objective function 1: Total Mass
    f1 = 3*(m_act + m_link);
    
    % Objective function 2: GCI
    gci = SPM_kinematics(diam_base/2,diam_plat/2,x(i,1),x(i,2),[],[]);
    if isValid
        f2 = -gci;
    else
        f2 = -(gci/penalty);
    end

    % Save the values to a table
    T(i,:) = [f1, f2, max(Torque), max(Force.Data(end,:)), m_act, m_link, deflection_mm];
    
    
    % fprintf('Outputs: m = %.2f, GCI = %.4f, T = %.2f, F = %.1f, m_act = %.1f, m_link = %.3f, defl = %.3f\n',...
    %     f1, f2, max(Torque), max(max(Force)), m_act, m_link, deflection_mm)
    if ~isValid
        fprintf('Penalty applied: %.1f\n',penalty)
    end
    catch ME
        warning("Simulation failed: %s", ME.message);
        clear mex;
    end

    if mod(i, 150) == 0
        bdclose('SPM_Model');    % Close the model
        load_system('SPM_Model'); % Re-open it clean
        set_param('SPM_Model', 'FastRestart', 'on'); % Re-enable fast restart
    end
end

% Save the results to a file
data = array2table(horzcat(x,T), 'VariableNames',{'h','theta', 'diam_plat','mass','gci','T','F','m_act','m_link','deflection_mm'});
writetable(data,'spm_halton.txt')

%% Visualize the results
pltSave = true;
close all;

% Open the saved solution
data = readtable('spm_halton.txt');

% Create a meshgrid for plotting
xq = linspace(min(data.h), max(data.h), 100);
yq = linspace(min(data.theta), max(data.theta), 100);
[X,Y] = meshgrid(xq, yq);
Z_mass = griddata(data.h, data.theta, data.mass, X, Y, 'cubic');
Z_gci = griddata(data.h, data.theta, data.gci, X, Y, 'cubic');

% Plot the design variables against mass
fg = figure(1);
mesh(X,Y,Z_mass)
hold on
plot3(data.h, data.theta, data.mass, '*')
hold off
xlabel('Height (m)')
ylabel('theta (deg)')
zlabel('Mass (g)')
title('SPM Mass Sensitivity')
if pltSave
    saveas(fg,'SPM_mass_sensitivity_halton.png')
end

% Plot the design variables against gci
fg = figure(2);
mesh(X,Y,Z_gci)
hold on
plot3(data.h, data.theta, data.gci, '*')
hold off
xlabel('Height (m)')
ylabel('theta (deg)')
zlabel('GCI')
title('SPM GCI Sensitivity')
if pltSave
    saveas(fg,'SPM_gci_sensitivity_halton.png')
end