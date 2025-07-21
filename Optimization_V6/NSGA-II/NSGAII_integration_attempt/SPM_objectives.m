function F = SPM_objectives(x, params)
% Enhanced SPM_objectives function with soft constraints and trade-off exploration
%
% Inputs:
%   x - Design vector [W, T, actuator_index]
%   params - Struct containing optimization parameters
%
% Outputs:
%   F - Objective vector [total_mass, -GCI]

% Deflection optimization with initial parameters
[mass, isValid, W, T] = Deflection_Optim(params.load, params.h, params.theta_deg, ...
                                         params.platform_diam, params.debug, x(1), x(2));

% Get actuator mass based on selected index
actuator_mass = params.actuator_masses(round(x(3)));

% Compute total mass
total_mass = mass + actuator_mass;

% Compute Global Conditioning Index (GCI)
% Use silent version to avoid excessive output
[gci] = SPM_kinematics_silent(params.base_rad, params.platform_diam/2, ...
                               params.h, params.theta_deg, [], [], isValid);

% Introduce a softer penalization strategy
if ~isValid
    % Smooth penalty function that doesn't completely reject designs
    % but gradually increases mass penalty for less valid designs
    gci_penalty = 1 - exp(-gci);  % Smooth penalty function
    total_mass = total_mass * (1 + gci_penalty);
end

% Return objectives: [total_mass, -GCI]
% Negative GCI because NSGA-II maximizes the second objective
F = [total_mass, -gci];
end