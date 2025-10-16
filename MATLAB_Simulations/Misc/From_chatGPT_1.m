% From_chatGPT_1

clear
clc

model = createpde();
geometryFromEdges(model, @circleg);

mesh = generateMesh(model);

[p, e, t] = meshToPet(mesh);

% Coordinates of points
x = p(1, :);
y = p(2, :);

% Find elements close to the center (radius threshold)
center_radius = 0.2; %0.2; % Adjust this for desired density
distances = sqrt(x.^2 + y.^2);
center_elements = distances < center_radius;

% Refine elements near the center
% mesh = refinemesh(model.Geometry, p, e, t, center_elements);
mesh = generateMesh(model, 'GeometricOrder', 'linear');

figure(11)
clf
pdeplot(model, 'Mesh', 'on');










