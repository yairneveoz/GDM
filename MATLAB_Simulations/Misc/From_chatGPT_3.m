% From_chatGPT_3

clear
clc

%%%
model = createpde();

% Define the outer and inner circles
R_outer = 1;   % Radius of the outer circle
R_inner = 0.2; % Radius of the inner circle

outer_circle = [1; 0; 0; R_outer; 0; 0; 0]; % Outer circle: [1; center_x; center_y; radius; 0; 0; 0]
inner_circle = [1; 0; 0; R_inner; 0; 0; 0]; % Inner circle: [1; center_x; center_y; radius; 0; 0; 0]

% Combine the circles into a geometry description matrix
gd = [outer_circle, inner_circle]; % Geometry description matrix

% Define the set formula (union of the two circles)
sf = 'C1+C2'; % Union: Keep both circles

% Define the labels for the geometry objects
ns = char('C1', 'C2')'; % Names for the geometry objects

% Create the composite geometry
geometry = decsg(gd, sf, ns);
geometryFromEdges(model, geometry);

% Visualize the geometry
figure(12);
pdegplot(model, 'FaceLabels', 'on'); % Show face labels to identify regions
title('Composite Geometry with Region Labels');

%%%
% Generate the mesh with region-specific Hmax
mesh = generateMesh(model, ...
    'Hmax', 0.1, ...  % Default maximum element size
    'Hedge', {(1), (2)});
% Create









