function [mass, actuator] = actuator_selection(torque, diam_plat)
% Given an actuator torque requirement, selects the lightest actuator that
% suffices. All actuators on the list already satisy necessary speed.
%
%   VARIABLES
%   ---------
%
%   torque := torque required from the central axis in N-m (scalar). This
%   should be the output value from the SPM model. We will compute from
%   this the necessary torque at the actuator.
%
%   diam_plat := platform diameter in m
%
%
%   OUTPUT
%   ---------
%
%   m := mass of the actuator in grams
%
%   actuator := table with the selected actuator and all parameters
%
%
%   USAGE EXAMPLE
%   -----
%
%   [m,~] = actuator_selection(200)
% 

% Import the motors spreadsheet and sort by mass
motors = readtable("motor_spreadsheet.csv");
motors = sortrows(motors,"Mass_g");

% Compute the required torque
gear_ratio = 6;
fos = 2;
max_t = max(torque)*fos/gear_ratio;


% Find the lightest actuator that meets the torque requirement
rowIdx = find(motors.MaxTorque_Nm > max_t, 1, 'first'); % Find first index

if isempty(rowIdx)
    disp('No motors have sufficient torque.');
    mass = max(motors.Mass_g)*3;
else
    actuator = motors(rowIdx, :);
    mass = actuator.Mass_g;
end
