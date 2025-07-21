function [GCI] = SPM_kinematics(base_rad,top_rad,h,theta,alpha1,alpha2,plotting)
% Given an h and theta OR alpha1 and alpha2, this function will output a conditioning
% index of a coaxial Spherical Parallel Manipulator.
%
%   VARIABLES
%   ---------
%
%   base_rad := radius of the base platform (scalar)
%
%   top_rad := radius of the mobile platform (scalar)
%
%   h := height of the SPM (scalar)
%
%   theta := angle made by the revolute joint axis and projection of
%   proximal-distal joint axis onto the mobile platform's plane (scalar)
% 
%   alpha1 := angle made between the base joint axis and the
%   proximal-distal joint axis (scalar) 
%
%   alpha2 := angle made between the revolute joint axis and the
%   proximal-distal joint axis (scalar)
%
%   plotting := boolean indicating if we are plotting the output.
%
%   USAGE
%   -----
% EX:
%   [GCI] = SPM_kinematics(base_rad,top_rad,h,theta,[],[])
%         
%   [GCI] = SPM_kinematics(base_rad,top_rad,[],theta,alpha1,[])
%
%   [GCI] = SPM_kinematics(base_rad,top_rad,[],theta,alpha1,alpha2)
%
% alpha2 is the only variable that is not necessary, the function requires
% a theta,base_rad, AND top_rad to run. h and alpha1 are interchangeable in
% necessity

%%

% renaming theta to a different angle so there's no confusion
unk = theta;

% phi, theta, and psi represent iterative rotations around X,Y,Z axes of
% the SPM
phi=0:15:360;
theta=0:15:360;
psi=0:15:360;

eta11 = 0;
eta12 = 120;
eta13 = 240;
gamma2 = 120;
beta2 = asind((2.*sqrt(3)./3) .* sind(gamma2./2));

if isempty(alpha1)
    alpha1 = atand(base_rad/h);
    % disp(alpha1)
elseif isempty(h)
    h = base_rad/(tand(alpha1));
else
    error('Need h or alpha1 to run')
end

if nargin < 7 || isempty(plotting)
    plotting = false;
end

%% Calculations for alpha2
if isempty(alpha2) && ~isempty(unk)
    a_mag = sqrt(h.^2 + base_rad.^2);
    
    Ry_clockwise = [cosd(alpha1), 0, sind(alpha1);
                    0,           1,           0;
                   -sind(alpha1),0,  cosd(alpha1)];
    
    Rz_counter = [cosd(unk), -sind(unk), 0;
                  sind(unk),  cosd(unk), 0;
                  0,          0,         1];
    
    a_vec = mtimes(Ry_clockwise, [0 0 -1]');
    a_vec = a_vec .* a_mag;
    
    x_vec = [1 0 0];
    
    b_vec = mtimes(Rz_counter,x_vec') .* top_rad; 
    
    alpha2 = acosd((dot(a_vec,b_vec)) ./ (norm(a_vec) .* norm(b_vec)));
    % disp(alpha2)
    if alpha2 < 0
        alpha2 = alpha2 + 360;
    end
end
%%
if isempty(unk)
    error('Need theta value')
end
%%
[W,N,nn,CN]=Workspace1(alpha1,eta11,eta12,eta13,alpha2,beta2,gamma2,phi,theta,psi,0.2);

% if ~isempty(W)
%     figure
%     plot3(W(1,:),W(2,:),W(3,:),'.')
%     xlabel("\phi")
%     ylabel("\theta")
%     zlabel("\psi")
%     %view(90,0)
%     grid on
%     grid minor
%     set(gcf,'color','w')
% end

%% Plot the results
if ~isempty(N)
    if plotting
        figure
        plot3(N(1,:),N(2,:),N(3,:),'.')
        hold on
        if ~isempty(nn)
            plot3(nn(1,:),nn(2,:),nn(3,:),'.')
        end
        hold on
        xlabel("x")
        ylabel("y")
        zlabel("z")
        zlim([0,1])
        grid on
        grid minor
        set(gcf,'color','w')
        x= N(1,:)';
        y= N(2,:)';
        z= N(3,:)';
        y = y(z>=0);
        x = x(z>=0);
        CN = CN(z>=0);
        z = z(z>=0);
        legend('non-singular','singular')
        str_title = sprintf(['Range of Motion: h = %d mm,',char(952),' = %d', char(176)],h,unk);
        title(str_title)
    
        figure
        
        A = -1:0.1:1;
        colorpost = CN';
        cmap = winter(length(CN));
        numcol = size(cmap,1);
        cps = floor(rescale(colorpost, 1, numcol));
        cps2 = cps./ max(cps);
        
        scatter3(x,y,z,5.*ones(size(x,1),1),colorpost,'filled')
        c=colorbar();
        ylabel(c,'Conditioning #','FontSize',8,'Rotation',270)
        xlim([-1 1])
        ylim([-1 1])
        zlim([0 1])
        grid on
        grid minor
        set(gcf,'color','w')
        xlabel("x")
        ylabel("y")
        zlabel("z")
        clim([0, 1]);
        str_title2 = sprintf(['Range of Motion w/ GCI: h = %d mm,',char(952),' = %d', char(176)],h,unk);
        title(str_title2)
    end
end

GCI = mean(CN);

end