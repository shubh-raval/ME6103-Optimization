function [T,F] = SPM_Simulation(h,theta,diam_base,diam_plat)
% Runs the Simulink SPM model to determine actuator torque and link force.
% Given an h, theta, and diameters, this function will output call the 
% Simulink multibody simulation and output the actuator torque and linkage 
% force as timeseries.
%
%   VARIABLES
%   ---------
%
%   h := height of the SPM in m (scalar)
%
%   theta := angle (degrees) made by the revolute joint axis and projection 
%   of proximal-distal joint axis onto the mobile platform's plane (scalar)
%
%   diam_base := base diameter in m
%
%   diam_plat := platform diameter in m
%   
%
%   OUTPUT
%   ---------
%
%   T := Peak Actuator Torque in N*m (timeseries)
%
%   F := Platform Linkage Force magnitude in N (timeseries). Each of the
%   data columns represents a different platform linkage.
%
%
%   USAGE
%   -----
% EX:
%   [T,F] = SPM_Simulation(.1, 90, .18, .14)
%
% 

%%

% Design Parameters
link_h = h*1000; %mm
link_theta = deg2rad(theta);
radius = diam_base*500; %mm

% Robot and Workpiece parameters
motorlink.width = 6; %mm
motorlink.length = radius+motorlink.width/2; %mm
motorlink.thickness = 6; %mm
motorlink.height = 8; %mm
motorlink.diam = 12; %mm
motorlink.angle = atan2(link_h, radius); %pi/2 - 

platform.diam = diam_plat*1000; %mm
platform.thickness = 8; %mm

% Pass the simulation parameters to the simulation
simIn = Simulink.SimulationInput('SPM_Model');
simIn = simIn.setVariable('motorlink', motorlink);
simIn = simIn.setVariable('platform', platform);
simIn = simIn.setVariable('link_h', link_h);
simIn = simIn.setVariable('link_theta', link_theta);
simIn = simIn.setVariable('radius', radius);

% Run the simulation and parse the outputs
simout = sim(simIn);
T = simout.yout{1}.Values;
F = simout.yout{2}.Values;