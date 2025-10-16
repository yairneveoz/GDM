% Color_density

clear
clc

size_i = 800; % pixels
size_j = 1000; % pixels

xx = 1:size_j;
yy = 1:size_i;

[XX,YY] = meshgrid(xx,yy);

color_0 = 0.5;
sigma = 75; % 180

% Set center point:
c1_j1 = ceil(size_j/4);
c1_i1 = ceil(size_i/2);
c1_j2 = ceil(size_j*3/4);
c1_i2 = ceil(size_i/2);

array01 = zeros(size_i,size_j);
array02 = zeros(size_i,size_j);
array01(c1_i1, c1_j1) = 1;
array02(c1_i2, c1_j2) = 1;

% array11 = imgaussfilt(array01, sigma, 'Padding', 'replicate');
% array12 = imgaussfilt(array02, sigma, 'Padding', 'replicate');
array11 = 1*exp(-0.5*(((XX - c1_j1).^2) + ((YY - c1_i1).^2))./sigma^2);
array12 = 1*exp(-0.5*(((XX - c1_j2).^2) + ((YY - c1_i2).^2))./sigma^2);

min_array1 = min(array11(:));
max_array1 = max(array11(:));

density_array_p1 = array11 + 1;
density_array_p2 = array12 + 1;

density_array_n2 = 1/1 - 1./density_array_p2;

density_array = density_array_p1 + density_array_n2 - 1;

norm_array1 = (array11 - min_array1)/(max_array1 - min_array1);
% norm_array2 = 1 - norm_array1;
norm_array2 = norm_array1;

[grad_x, grad_y] = gradient(norm_array2);

figure(7)
clf
% subplot(3,1,1:2)
imagesc(density_array_p1)
hold on
contour(density_array_p1,'r')
hold off
colormap(gray)
colorbar
caxis([0 2])
axis equal
axis tight

