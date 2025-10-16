% Universe_space_density_pixels

clear
clc
%%
array_size = 300;
rho_array0 = 0.5*ones(array_size);
weight0 = array_size^2*0.5;
N = 1000;
%
%% Mask
mask_size = 15;
normalized_mask = getMask(mask_size);
%
%%
linind0 = randsample(array_size^2, N, 'false');
point_array0 = zeros(array_size);
point_array0(linind0) = 1;
point_array = point_array0;
rho_array = imfilter(point_array, normalized_mask, 'replicate');
rho_array = rho_array + rho_array0;
%
%% 
weight = sum(rho_array(:));
rho_array = (rho_array - 1)*(weight/weight0) + 1;
disp(min(rho_array(:)));
figure(1)
clf
image(255*rho_array)
axis equal
axis tight
colormap(1-gray)
colorbar
%
%% Get gradient of density:
figure(2)
[point_array_x, point_array_y] = find(point_array);
[grad_rho_array_x,grad_rho_array_y] = gradient(rho_array);

% grad_point_array_x = zeros(array_size);
% grad_point_array_y = zeros(array_size);
grad_point_array_x = grad_rho_array_x(linind0);
grad_point_array_y = grad_rho_array_y(linind0);

quiver(point_array_x, point_array_y,...
    grad_point_array_x,grad_point_array_y)
axis equal
axis tight

figure(3)
grad_magnitude = sqrt(grad_rho_array_x.^2 + grad_rho_array_y.^2);
imagesc(grad_magnitude)
axis equal
axis tight