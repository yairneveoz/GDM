% Two_dimensional_array_of_springs

%{
Tracks the dynamics (locations over time) of an array of
points connected by 'springs' in one dimension.
%}

clear
clc

blue_white_red = flipud(blueWhiteRedColormap(64));

dx = 1;
dy = dx;

%{
Distance between two points when the spring is at 
its resting length
%}
dt = 0.01; % sec.
Nx = 31; % Number of points in x direction.
Ny = Nx;

c_T = 3; % Transverse light velocity.
c_L = c_T*1.608; % Longitudinal light velocity.


x0 = -(Nx/2)*dx:dx:(Nx/2)*dx;
y0 = -(Ny/2)*dx:dx:(Ny/2)*dy;

[X0, Y0] = meshgrid(x0,y0);

X = X0;
Y = Y0;

ux = 0.3;
uy = 0.3;
rho_0 = 1;

center_x = ceil(Nx/2);
center_y = ceil(Ny/2);

% X(center_x,center_y) = X(center_x,center_y) + ux;
% Y(center_x,center_y) = Y(center_x,center_y) + uy;
% Y(center_x,center_y+1) = Y(center_x,center_y+1) + 1.6*uy;

% Y(:,center_y) = Y0(:,center_y) + uy*sin(2*Y0(:,center_y));
sigma_x = 2;
sigma_y = sigma_x;

R0 = sqrt(X0.^2 + Y0.^2);

% X = X0 + ux*sin(-0.5*X0);
% Y = Y0 + uy*sin(-0.7*Y0);
X = X0 - sign(X0)*20*ux./(1+R0);
Y = Y0 - sign(Y0)*20*uy./(1+R0);


[V, C] = voronoin([X(:), Y(:)]);
%%%
figure(3)
clf
hold on
if 1
    % Calculate area of each Voronoi cell
    areas = zeros(length(C), 1); % Initialize an array for areas
    for i = 1:length(C)
        % Get the vertices of the ith Voronoi cell
        cellVertices = V(C{i}, :);

        % Check if the cell is bounded
        if all(C{i} ~= 1)  % Skip unbounded cells (infinite)
            areas(i) = polyarea(cellVertices(:, 1), cellVertices(:, 2));
        else
            areas(i) = NaN; % Mark unbounded cells with NaN
        end
        patch(V(C{i},1), V(C{i},2), areas(i))
    end
end
% areas(isnan(areas)) = 0;
scatter(X(:), Y(:), 4, 'k', 'filled')
colormap(blue_white_red)
axis([-Nx/2 Nx/2 -Nx/2 Nx/2])
caxis([0 2])
hold off
axis equal
% axis tight
%
%%
%{
% Calculate polygon areas
areas = polyarea(VoronoiDiagram);

% Define a colormap (e.g., jet)
cmap = colormap(jet);

% Interpolate colormap values based on polygon areas
num_colors = length(cmap);
color_indices = ceil(areas * (num_colors - 1));

% Assign colors to polygons
colors = cmap(color_indices, :);

% patch(V(C{1,:}), C, areas)
% figure(4)
% voronoi(X(:,1),X(:,2))
% voronoi(X(:), Y(:));
% hold on

% for i = 1:size(X(:),1)
%     A=V(C{i},:); %A=V(R{i},:);
%     B=A(any(~isinf(A),2),:);
%     if(size(B,1)>2)
% %         plot(polyshape(B),'FaceColor','b');
%         patch(B(:,1),B(:,2),'b');
%     end
% end
%}
