% Frame_dragging_with_voronoi_cells
%{
Create mesh.
Import mesh points.
Draw Voronoi cells.
Define center.
Rotate around center with steps over time.
Redraw Voronoi cells at every step.
Calculate bending of the lines.
Look for limit.
Do the same dilated and compressed space.
%}
clc
clear

locations = csvread('initial_locations1.csv');
x0 = 10*locations(:,1);
y0 = 10*locations(:,2);

r_min = 50;
r_max = 240;
[theta0, r0] = cart2pol(x0,y0);
theta0_r0 = [theta0, r0];
in_r_min = theta0_r0(:,2) < r_min;
theta0_r0(in_r_min,:) = [];
theta1 = theta0_r0(:,1); 
r1 = theta0_r0(:,2); 

phi = linspace(0,2*pi,20);
circ_x = r_min*cos(phi);
circ_y = r_min*sin(phi);
%
%%
L = 20; % L = m*v cross r, v = r*d(theta)/d(t)
delta_t = 0.1;
delta_theta = L./(r1);
theta = theta1 + delta_theta;

[x, y] = pol2cart(theta,r1);

figure(4)
theta = theta1;
for iter = 1:200
    delta_theta = delta_t*(L./r1 - L./r_max);
    theta = theta + delta_theta;
    
    [x, y] = pol2cart(theta,r1);
%     [vx, vy] = voronoi(x, y);

    % Plot Voronoi diagram
    subplot(1,2,1)
    plot(circ_x, circ_y, '-')
    hold on
    voronoi(x,y)
    hold off
    axis equal
    axis(r_max*[-1 1 -1 1])
    title(num2str(iter))
    pause(0.05)
    drawnow
    
%         subplot(1,2,2)
%     plot(circ_x, circ_y, '-')
%     hold on
%     plot(vx, vy, 'b-');
%     hold off
%     axis equal
    
    axis(r_max*[-1 1 -1 1])
    title(num2str(iter))
    pause(0.1)
    drawnow
end


% voronoi(x,y)
% hold off
% axis equal
%
%%



