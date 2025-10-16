% function [spiral_path] = spirals(r,R,v)
clear
clc

n_r = 10;
n_R = 11;
omega_r = 50;
omega_R = 1;
t = 0:0.01:10;
phi_r = omega_r*t;

r = 1; % /gamma
R = 20; %137.036; % /gamma
v_r = 0.8;
v_R = 5;

e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

x_R = R*cos(omega_R*t); % /gamma
y_R = R*sin(omega_R*t); % /gamma
z_R = v_R*t; % /gamma

R_t = [x_R', y_R', z_R'];
radial_R = [x_R',y_R',zeros(size(x_R'))];

diff_x_R = diff(x_R);
diff_y_R = diff(y_R);
diff_z_R = diff(z_R);
diff_R = [diff_x_R', diff_y_R', diff_z_R'];

%%%
x_r = r*cos(omega_r*t);
y_r = r*sin(omega_r*t);
z_r = 0*v_r*t;
%%%

r_t = [x_r', y_r', z_r'];

x_rR = x_r + x_R;
y_rR = y_r + y_R;
z_rR = z_r + z_R;

rR_t = R_t + r_t;

nt = 100;
s = R;
if 1
    
    % Initial:
    ni = 1;
    start_e1 = diff_R(ni,:)/norm(diff_R(ni,:));
    start_e2 = radial_R(ni,:)/norm(radial_R(ni,:));
    start_e3 = cross(start_e2,start_e1);
    
    % After rotation:
    rot_e1 = diff_R(nt,:)/norm(diff_R(nt,:));
    rot_e2 = radial_R(nt,:)/norm(radial_R(nt,:));
    rot_e3 = cross(rot_e2,rot_e1);

    figure(40)
    clf
    subplot(2,2,1)
    plot3(R_t(:,1), R_t(:,2), R_t(:,3), '-')

%     plot3(x_R, y_R, z_R, '-')
    hold on
    % Tangent:
    plot3([R_t(1,1), R_t(1,1) + s*start_e1(1)],...
        [R_t(1,2), R_t(1,2) + s*start_e1(2)],...
        [R_t(1,3), R_t(1,3) + s*start_e1(3)],'r-','LineWidth',1)
    
    % Radial:
    plot3([R_t(1,1), R_t(1,1) + s*start_e2(1)],...
        [R_t(1,2), R_t(1,2) + s*start_e2(2)],...
        [R_t(1,3), R_t(1,3) + s*start_e2(3)],'g-','LineWidth',1)
    
    % Orthogonal to:
    plot3([R_t(1,1), R_t(1,1) + s*start_e3(1)],...
        [R_t(1,2), R_t(1,2) + s*start_e3(2)],...
        [R_t(1,3), R_t(1,3) + s*start_e3(3)],'b-','LineWidth',1)
    
    % R
    plot3([0 R_t(nt,1)],[0 R_t(nt,2)],[R_t(nt,3) R_t(nt,3)],'k-')

    % Tangent:
    plot3([R_t(nt,1), R_t(nt,1) + s*rot_e1(1)],...
        [R_t(nt,2), R_t(nt,2) + s*rot_e1(2)],...
        [R_t(nt,3), R_t(nt,3) + s*rot_e1(3)],'r-','LineWidth',2)

    % Radial:
    plot3([R_t(nt,1), R_t(nt,1) + s*rot_e2(1)],...
        [R_t(nt,2), R_t(nt,2) + s*rot_e2(2)],...
        [R_t(nt,3), R_t(nt,3) + s*rot_e2(3)],'g-','LineWidth',2)

    % Orthogonal to:
    plot3([R_t(nt,1), R_t(nt,1) + s*rot_e3(1)],...
        [R_t(nt,2), R_t(nt,2) + s*rot_e3(2)],...
        [R_t(nt,3), R_t(nt,3) + s*rot_e3(3)],'b-','LineWidth',2)
    

    hold off
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
    
    subplot(2,2,2)
    plot3(x_r, y_r, z_r, 'r-')
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
    
    
    subplot(2,2,3)
    plot3(x_R, y_R, z_R, 'b-')
    hold on
    plot3(x_rR, y_rR, z_rR, 'r-')
    hold off
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
    
    subplot(2,2,4)
    plot3([0 rot_e1(1)],[0, rot_e1(2)], [0, rot_e1(3)],'r-')
    hold on
    plot3([0 rot_e2(1)],[0, rot_e2(2)], [0, rot_e2(3)],'g-')
	plot3([0 rot_e3(1)],[0, rot_e3(2)], [0, rot_e3(3)],'b-')
    hold off
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
end
% end