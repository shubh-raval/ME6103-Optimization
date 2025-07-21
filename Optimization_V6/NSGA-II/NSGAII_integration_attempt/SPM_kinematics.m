function [GCI, x_data, y_data, z_data, CN] = SPM_kinematics(base_rad, top_rad, h, theta, alpha1, alpha2, isValid)
% SPM_kinematics computes the Global Conditioning Index (GCI) for a coaxial
% Spherical Parallel Manipulator (SPM) given either h or alpha1 (or both).
%
%   [GCI, x_data, y_data, z_data, CN] = SPM_kinematics(base_rad, top_rad, h, theta, alpha1, alpha2, isValid)
%
%   VARIABLES
%   ---------
%   base_rad := radius of the base platform (scalar)
%   top_rad  := radius of the mobile platform (scalar)
%   h        := height of the SPM (scalar)
%   theta    := angle made by the revolute joint axis and projection of the 
%               proximal-distal joint axis onto the mobile platform's plane (scalar)
%   alpha1   := angle made between the base joint axis and the proximal-distal joint axis (scalar)
%   alpha2   := angle made between the revolute joint axis and the proximal-distal joint axis (scalar)
%   isValid  := boolean indicating if a valid linkage was found. If false,
%               a penalty factor is applied to the calculated GCI.
%
%   USAGE EXAMPLES:
%       [GCI, x_data, y_data, z_data, CN] = SPM_kinematics(base_rad, top_rad, h, theta, [], []);
%       [GCI, x_data, y_data, z_data, CN] = SPM_kinematics(base_rad, top_rad, [], theta, alpha1, []);
%       [GCI, x_data, y_data, z_data, CN] = SPM_kinematics(base_rad, top_rad, [], theta, alpha1, alpha2);
%
%   Note: h and alpha1 are interchangeable. If both are provided,
%         the function will use the supplied values.

% Check for required inputs and compute missing ones if possible
if isempty(alpha1) && ~isempty(h)
    alpha1 = atand(base_rad/h);
elseif isempty(h) && ~isempty(alpha1)
    h = base_rad/(tand(alpha1));
elseif isempty(alpha1) && isempty(h)
    error('Need h or alpha1 to run');
end

% Display input values for debugging
fprintf('SPM_kinematics called with h = %f, alpha1 = %f\n', h, alpha1);

% Rename theta to unk to avoid confusion
unk = theta;

% Define iterative angles for workspace sampling (in degrees) using a coarser grid:
phi = 0:15:360;
theta = 0:15:360;
psi = 0:15:360;

% Base joint angular offsets (degrees)
eta11 = 0;
eta12 = 120;
eta13 = 240;

% Fixed parameters for platform kinematics
beta2 = 90;
gamma2 = 120;

%% Calculations for alpha2 (if not provided)
if isempty(alpha2) && ~isempty(unk)
    a_mag = sqrt(h.^2 + base_rad.^2);
    
    Ry_clockwise = [cosd(alpha1), 0, sind(alpha1);
                    0,           1,           0;
                   -sind(alpha1),0,  cosd(alpha1)];
    
    Rz_counter = [cosd(unk), -sind(unk), 0;
                  sind(unk),  cosd(unk), 0;
                  0,          0,         1];
    
    a_vec = Ry_clockwise * [0; 0; -1];
    a_vec = a_vec * a_mag;
    
    x_vec = [1; 0; 0];
    
    b_vec = Rz_counter * x_vec * top_rad;
    
    alpha2 = rad2deg(dot(a_vec, b_vec) / (norm(a_vec) * norm(b_vec)));
    if alpha2 < 0
        alpha2 = alpha2 + 360;
    end
end

%% Check for required theta value
if isempty(unk)
    error('Need theta value');
end

%% Evaluate Workspace and Conditioning Index (GCI)
[W, N, nn, CN] = Workspace1(alpha1, eta11, eta12, eta13, alpha2, beta2, gamma2, phi, theta, psi, 0.2);

% Extract data for plotting
x_data = N(1,:)';
y_data = N(2,:)';
z_data = N(3,:)';

idx = z_data >= 0;
x_data = x_data(idx);
y_data = y_data(idx);
z_data = z_data(idx);
CN = CN(idx);

% Calculate mean Global Conditioning Index (GCI)
GCI = mean(CN);

% Apply penalty factor for invalid designs
if ~isValid
   GCI = GCI * 0.1;
end
end