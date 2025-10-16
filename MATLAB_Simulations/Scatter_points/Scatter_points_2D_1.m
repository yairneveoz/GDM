% Scatter_points_2D_1

clear
clc
close all

% set(0,'DefaultFigureRenderer','painters');   % safest for 2D

% Initialize random seed for reproducibility
% rng('shuffle');
N = 200;
L = 1; % Define the limits for the scatter plot
a = 0; % Define the lower limit for the y-axis
x = L*rand(N, 1); % Generate N random x-coordinates
y = L*rand(N, 1); % Generate N random y-coordinates

l0 = sqrt(L^2/N);

figure('Color','w','Position',[100 100 640 480],'Visible','on');
for iter = 1:1
    %% Generate Voronoi
    % Compute Delaunay triangulation
    % [vx, vy] = voronoi(x, y);
    [v, c] = voronoin([x, y]);
    

    new_x = zeros(size(x));
    new_y = zeros(size(y));
    
    % for ind = 1:length(c)
    %     vx_ind = v(c{ind}, 1);
    %     vy_ind = v(c{ind}, 2);
    % 
    %     vx_center = mean(vx_ind);
    %     vy_center = mean(vy_ind);
    % 
    %     new_x(ind) = vx_center;
    %     new_y(ind) = vy_center;
    % end
    
    
    voronoi(x, y);
    hold on    
    % voronoi(new_x, new_y);
    hold off
    xlabel('X');
    ylabel('Y');
    title('Voronoi Diagram of Random Points');
    axis([0 L a L])
    axis equal
    % axis tight
    box on
end
