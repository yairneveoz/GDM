% random_points_dynamics_2D

clear
clc

N = 200; % number of points
L = 10; % smaple length

x0 = L*rand(N,1);
y0 = L*rand(N,1);

x = x0;
y = y0;

L_factor = 1;
T_factor = 1;
dx0 = L/N;
dy0 = L/N;

% DX0 = dx0*ones(size(x));
DX0 = dx0*0.5*randn(size(x));
DY0 = dy0*0.5*randn(size(x));

figure(4)
subplot(2,2,1)
scatter(x, y, '.')
axis([0 L 0 L])
axis equal
axis tight
xlabel('X')
ylabel('Y')

%%%
[v,c] = voronoin([x, y]);
subplot(2,2,2)
voronoi(x, y)
axis([0 L 0 L])
axis equal
axis tight
xlabel('X')
ylabel('Y')

%% %

xmin = 0; xmax = L;
ymin = 0; ymax = L;
box = polyshape([xmin xmax xmax xmin], [ymin ymin ymax ymax]);

% Generate Voronoi diagram
points = [x, y];
[V, C] = voronoin(points);

% Loop through cells
for i = 1:length(C)
    verts = V(C{i}, :);
    if all(isfinite(verts(:)))  % Already bounded
        pgon = polyshape(verts(:,1), verts(:,2));
    else
        % Infinite cell: try clipping
        verts = remove_infinite_edges_and_extend(verts, box); % custom function
        pgon = intersect(polyshape(verts), box);  % Clip with box
    end
    
    area = area(pgon);
    % Store or use area
end

%%
x = cellVertices(:, 1);
y = cellVertices(:, 2);
area = polyarea(x, y);

%%
subplot(2,2,3)
voronoi(x, y)
axis([0 L 0 L])
axis equal
axis tight
xlabel('X')
ylabel('Y')

%%
for iter = 1:1000
    disp(iter)
    %%%
    mean_diff_x = zeros(size(x));
    for n = 2:N-1
        mean_diff_x(n) = mean([abs(x(n) - x(n-1)),...
            abs(x(n) - x(n+1))]);
    end
    mean_diff_x(1) = abs(x(2) - x(1));
    mean_diff_x(N) = abs(x(N) - x(N-1));
    %%%

    
    
    xlabel('x')
    title('Initial spatial distribution')


    subplot(3,1,2)
    plot(mean_diff_x, '.-')
%     plot(mean_diff_x/dx0, '.-')
    % ylim([0 2])
    xlabel('Cell index')
    ylabel('Cell size')
    title('Cell sizes by index')
    
    DX = DX0.*mean_diff_x;
%     x1 = x + DX;
    %%% Periodic
    x1_0 = mod(x + DX, L);
    x1 = sort(x1_0);
    %%%

    subplot(3,1,3)
    plot(x1, 1*ones(size(x0)), 'r.')
    ylim([0 2])
    pause(0.05)
    drawnow
    
    x = x1;
    
end