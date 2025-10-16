% Generate 200 random 2D points

clear
clc

rng(0);
L = 10;
points = L*rand(200, 2);

% Compute the Voronoi diagram
[V, C] = voronoin(points);

% Define the bounding box (unit square)
boxShape = polyshape(L*[0 1 1 0], L*[0 0 1 1]);

% Initialize area array
areas = nan(length(C), 1);

for i = 1:length(C)
    region = C{i};
    
    % Skip cells with infinite vertices
    if any(region == 1) || any(region == 0) || any(isinf(region)) || isempty(region)
        continue;
    end

    % Get the polygon vertices
    try
        cellVerts = V(region, :);
        if any(isinf(cellVerts(:)))
            continue;
        end

        % Create the cell polygon
        cellPoly = polyshape(cellVerts(:,1), cellVerts(:,2));
        
        % Clip the polygon to the bounding box
        clipped = intersect(cellPoly, boxShape);

        % Store the area
        areas(i) = area(clipped);

    catch
        continue;
    end
end

% Show results in table
T = table(points(:,1), points(:,2), areas, ...
    'VariableNames', {'X', 'Y', 'VoronoiCellArea'});
disp(T);

% Optional: Plot the diagram and cell areas
figure(4);
clf
plot(boxShape, 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'k');
hold on;
for i = 1:length(C)
    if isnan(areas(i))
        continue;
    end
    cellVerts = V(C{i}, :);
    pgon = polyshape(cellVerts(:,1), cellVerts(:,2));
    clipped = intersect(pgon, boxShape);
    plot(clipped, 'FaceColor', [0.5 0.8 1], 'EdgeColor', 'b');
end
plot(points(:,1), points(:,2), 'k.');
hold on
voronoi(points(:,1), points(:,2))
hold off
title('Voronoi Diagram (Clipped) with Colored Cells');
axis equal tight;
