% Random_Gaussian_points

clear
clc

N = 10000;
array_size_x = 600;
array_size_y = 500;

center_x = ceil(array_size_x/2);
center_y = ceil(array_size_y/2);

sigma = 80;

% rand_x = array_size_x*rand(N,1);
% rand_y = array_size_y*rand(N,1);
% rand_Pr = 1*rand(N,1);

randn_x = sigma*randn(N,1) + center_x;
randn_y = sigma*randn(N,1) + center_y;
rand_Pr = 1*rand(N,1);

figure(9)
plot(randn_x, randn_y, '.')
axis equal








