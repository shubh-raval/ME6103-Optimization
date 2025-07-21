clear all; 
close all;
clc; 
%% ME6103 Assignment 3 Shubh Raval


%% Design Variables
D_p = 0.5; %Pipe Diameter in m
L_e = 6; % Length of the enclosure
N = 3; % Number of turns starts at 0

%% Visualization
plotParallelLines(N,L_e,D_p) 

%% Dependent Variables
S_p = D_p; %Pipe Spacing = Pipe Diameter
H_e = D_p*1.1; %Slightly Larger than Pipe Diameter
W_e = N*(D_p+S_p) + D_p; % Width is the number of turns + 1 addition exit Pipe
Half_Circ = D_p*pi*.5;

%% Volume Calculation
Volume = W_e * H_e * L_e;
fprintf('The volume is %s \n', Volume)


%% Total Pipe Length
%Straight Run on the bottom Pipe and the Exit Pipe
ExposedStraightPipe_Len = (N+1)*((L_e - 1.5*D_p));  % There will be N+1 of these
if N>1
    ExposedStraightPipe_Len = ExposedStraightPipe_Len + (N-1)*(L_e-(3*D_p));
end

U_bend_Len = N*((1.5*D_p*pi*.5)); 


Total_Len =   U_bend_Len + ExposedStraightPipe_Len; % m
Total_surface_area = Total_Len * Half_Circ;

%% Required Heat Transfer Params and Constants
m_dot = .2 ; % kg/s
cp_water = 4184; % J*kg/K
q_sun = 1000; %W/m^2
rho_water = 998;

%% Assume pipe has 100% absortivity from sun so that removes conduction and thickness << D_p
Delta_T = 3; % Required minimum temp increase
Q = m_dot*cp_water*Delta_T;
minimum_SA = Q/q_sun;

fprintf('The minimum SA is %s \n',minimum_SA);
fprintf('The current SA is %s \n', Total_surface_area);

%% Total Headloss
f_friction = 0.04; 
flow_v = m_dot/(rho_water*pi*(D_p*.5)^2);
mhl_t1 = ((flow_v)^2)*(1/(2*9.81));
major_head_loss = f_friction* (ExposedStraightPipe_Len/D_p)*mhl_t1;

% Using Equivalent Piping Length 
ratio = 50; %Return Bend

k = f_friction* ratio;

% minor loss
minor_head_loss = k * ((flow_v)^2)/(2*9.81); %One Bend

total_head_loss = major_head_loss + N * minor_head_loss;

fprintf('The total Head loss is %s \n',total_head_loss);




%nsga_2(200,100)
% Read the solution file
data = dlmread('solution.txt');

% Extract f1 (column 4) and f2 (column 5)
f1_values = data(:, 4);
f2_values = data(:, 5);

% Compute max values
f1_max = max(f1_values);
f2_max = max(f2_values);
f1_min = min(f1_values);
f2_min = min(f2_values);
% Display the results
fprintf('f1_max = %.6f\n', f1_max);
fprintf('f2_max = %.6f\n', f2_max);
fprintf('f1_min = %.6f\n', f1_min);
fprintf('f2_min = %.6f\n', f2_min);




