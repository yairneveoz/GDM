% From_chatGPT_2


model = createpde();

% Define the outer and inner circles
R_outer = 1;   % Radius of the outer circle
R_inner = 0.2; % Radius of the inner circle

outer_circle = [1; 0; 0; R_outer]; % Circle: [1; center_x; center_y; radius]
inner_circle = [1; 0; 0; R_inner];
circles = [outer_circle, inner_circle];

% Create a set formula for the geometry
set_formula = '1+2'; % Union of two circles
geometry = decsg(circles, set_formula, char('C1', 'C2')');
geometryFromEdges(model, geometry);

% Generate a mesh with regional control
generateMesh(model, ...
    'Hmax', 0.1, ...           % Default maximum element size
    'Hedge', [0.02, 0.1]);    % Specify H for individual regions

figure(12);

% Visualize the mesh
pdeplot(model, 'Mesh', 'on');
title('Composite Mesh with Smaller Circle and Finer Mesh');

   





