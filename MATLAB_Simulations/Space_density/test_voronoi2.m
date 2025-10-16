% Generate initial random points

clear
clc

rng(0);
points = rand(200, 2);

% Compute Voronoi diagram
[V, C] = voronoin(points);
boxShape = polyshape([0 1 1 0], [0 0 1 1]);

% Precompute clipped Voronoi polygons for each point
cellPolys = cell(length(C), 1);

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
        end
    catch
        continue;
    end
end

% Animate movement inside cells
figure;
axis equal tight;
xlim([0 1]); ylim([0 1]);
hold on;

for t = 1:200  % number of animation steps
    clf;
    plot(boxShape, 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'k');
    hold on;
    
    for i = 1:length(points)
        % Try to move point randomly
        if isempty(cellPolys{i})
            continue;
        end
        move = 0.005 * randn(1, 2);  % small random step
        newPoint = points(i, :) + move;
        % Check if new point is inside its Voronoi cell
        if isinterior(cellPolys{i}, newPoint)
            points(i, :) = newPoint;
        end
        plot(points(i,1), points(i,2), 'k.', 'MarkerSize', 8);
    end
    
    title(['Voronoi Movement Step ', num2str(t)]);
    drawnow;
end
