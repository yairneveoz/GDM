    % Generate initial random points

clear
clc

rng(0);
L = 10;
points = L*rand(500, 2);

% Compute Voronoi diagram
[V, C] = voronoin(points);
boxShape = polyshape(L*[0 1 1 0], L*[0 0 1 1]);

% Precompute clipped Voronoi polygons and their areas
cellPolys = cell(length(C), 1);
cellStepSize = 1*nan(length(C), 1);

for i = 1:length(C)
    region = C{i};
    if isempty(region) || any(region == 1) || any(region == 0)
        continue;
    end
    try
        verts = V(region, :);
        if any(isinf(verts(:)))
            continue;
        end
        pgon = polyshape(verts(:,1), verts(:,2));
        clipped = intersect(pgon, boxShape);
        if clipped.NumRegions == 1
            cellPolys{i} = clipped;
            cellStepSize(i) = sqrt(area(clipped));  % step size ∝ sqrt(area)
        end
    catch
        continue;
    end
end

% Normalize step sizes for consistent animation scale
% cellStepSize = 0.01 * cellStepSize/max(cellStepSize);  % scale max step to ~0.01
% cellStepSize = 0.1 * cellStepSize/max(cellStepSize);  % scale max step to ~0.01
cellStepSize = 0.5*cellStepSize/L;  % scale max step to ~0.01

% Animate movement inside cells
figure(5);
axis equal tight;
xlim(L*[0 1]); ylim(L*[0 1]);
hold on;

for t = 1:1000 %200  % number of animation steps
    clf;
    plot(boxShape, 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'k');
    hold on;
    
    for i = 1:length(points)
        % Try to move point randomly within its cell
        if isempty(cellPolys{i}) || isnan(cellStepSize(i))
            continue;
        end
        move = cellStepSize(i) * randn(1, 2);  % random step scaled to cell size
        newPoint = points(i, :) + move;
        
        % Accept only if new point is inside its Voronoi cell
        if isinterior(cellPolys{i}, newPoint)
            points(i, :) = newPoint;
        end
        
        plot(points(i,1), points(i,2), 'k.', 'MarkerSize', 8);
    end
    voronoi(points(:,1), points(:,2))
    title(['Voronoi Cell Movement - Step ', num2str(t)]);
    drawnow;
end
