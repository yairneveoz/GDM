% Color_density_velocity

clear
clc

size_i = 800; % pixels
size_j = 1000; % pixels

color_0 = 0.5;
sigma = 250; % 180

% Set center point:
c1_j = ceil(size_j/2);
c1_i = ceil(size_i/2);

array0 = zeros(size_i,size_j);
array0(c1_i, c1_j) = 1;

array1 = imgaussfilt(array0, sigma, 'Padding', 'replicate');

min_array1 = min(array1(:));
max_array1 = max(array1(:));

norm_array1 = (array1 - min_array1)/(max_array1 - min_array1);
% norm_array2 = 1 - norm_array1;
norm_array2 = norm_array1;

[grad_x, grad_y] = gradient(norm_array2);
% 
%%

j_front_t0 = 100:1:900;
i_front_t0 = 1*ones(size(j_front_t0));

j_front_t = j_front_t0;
i_front_t = i_front_t0;

% Plot
figure(7)
clf
% subplot(3,1,1:2)
imagesc(norm_array2)
hold on
colormap(gray)
colorbar
caxis([0 1])
axis equal
axis tight

delta_t = 1; % 1.0

front_mask_linind = ...
    sub2ind([size_i, size_j], i_front_t0, j_front_t0); 

speed_factor = 1;
speed_data = speed_factor*norm_array2(front_mask_linind);

V_x = speed_factor*grad_x(front_mask_linind);
V_y = speed_factor*grad_y(front_mask_linind);
V_factor = 150;
%V_factor
%%
for t = 1:1000 % 700
%     
%     i_front_t = 1*(i_front_t + speed_data*delta_t);
%     j_front_t = 1*(j_front_t + speed_data*delta_t);

%     i_front_t = 1*(i_front_t + (speed_data - V_factor*V_y)*delta_t);
%     j_front_t = 1*(j_front_t + (-V_factor*V_x)*delta_t);
    i_front_t = 1*(i_front_t + (speed_data + V_factor*V_y)*delta_t);
    j_front_t = 1*(j_front_t + (+V_factor*V_x)*delta_t);


    hold on
    h1 = plot(j_front_t, i_front_t, 'g-');
    hold off
    pause(0.01)
    drawnow
    delete(h1)
    
    front_mask_linind = ...
        sub2ind([size_i, size_j], ...
        round(i_front_t), round(j_front_t));
    
    speed_data = speed_factor*norm_array2(front_mask_linind);
    
    V_x = speed_factor*grad_x(front_mask_linind);
    V_y = speed_factor*grad_y(front_mask_linind);
end

%
%%
%{
% [FX,FY] = gradient(F)
[grad_norm_array2_X, grad_norm_array2_Y] = ...
    gradient(norm_array2);

figure(8)
subplot(1,2,1)
imagesc(grad_norm_array2_X(:,200:800))
colorbar
subplot(1,2,2)
imagesc(grad_norm_array2_Y(100:700,:))
colorbar
%}


