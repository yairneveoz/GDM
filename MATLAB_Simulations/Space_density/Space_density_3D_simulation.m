% Space_density_3D_simulation
clear
clc

array_size = 200;
N = 2000;

X0 = array_size*rand(N,1);
Y0 = array_size*rand(N,1);
Z0 = array_size*rand(N,1);

% f1 = 5;
% f2 = -2;
% f3 = 0.001;
% r1 = 3;
% r2 = 10;
r1 = 3;
r2 = 15;
rr2 = 0.5:0.1:100;
%
%% plot
X = X0;
Y = Y0;
Z = Z0;

grad_scale = 0.1;
N_iterations = 100;

steps1 = [5,-2,1e-3,3,10];
steps2 = [2,-2,1e-3,3,10];
steps3 = [5,-2,1e-2,3,10];
steps4 = [5,-2,1e-4,3,10];
steps5 = [5,-2,1e-1,3,10];
steps6 = [5,-1,1e-3,3,10];
steps7 = [5,-2,0,3,15];
steps9 = [5,-2,0,3,15];
steps10 = [5,-1,0,3,15];
% for 3D:
steps11 = [5,-2,1e-3,3,15];
steps12 = [5,-1,1e-3,3,15];
s = steps11;

figure(2)
subplot(1,3,1)
FR = FRsteps(s(1),s(2),s(3),s(4),s(5),rr2);
plot(rr2, FR, '-','LineWidth',2)
hold on
plot([rr2(1), rr2(end)], [0 0], 'k--')
hold off
axis square
axis([0 20 -10 10])
xlabel('r')
ylabel('F(r)')
title('F(r)')
text(1,s(1)+1,num2str(s(1)))
text(s(4)+1,s(2)-1,num2str(s(2)))
text(s(5)+1,s(3)+1,num2str(s(3)))

for iter = 1:N_iterations
    dist_X = distD(X);
    dist_Y = distD(Y);
    dist_Z = distD(Z);
    dist_R = sqrt(dist_X.^2 + dist_Y.^2 + dist_Z.^2);
    if iter == 1
        subplot(1,3,2)
        plot3(X, Y, Z, '.')
        axis equal
        axis([0 array_size 0 array_size 0 array_size])
        title(['Iter = ',num2str(iter)])
        grid on
        drawnow
    end
    F_array = FRsteps(s(1),s(2),s(3),s(4),s(5),dist_R);

    deltaX_array = F_array.*dist_X./dist_R;
    deltaY_array = F_array.*dist_Y./dist_R;
    deltaZ_array = F_array.*dist_Z./dist_R;
    
    deltaX_array(isnan(deltaX_array)) = 0;
    deltaY_array(isnan(deltaY_array)) = 0;
    deltaZ_array(isnan(deltaZ_array)) = 0;
    
    deltaX = sum(deltaX_array,2)*grad_scale;
    deltaY = sum(deltaY_array,2)*grad_scale;
    deltaZ = sum(deltaZ_array,2)*grad_scale;

    subplot(1,3,3)
    plot3(X, Y, Z, '.')
    axis equal
    axis([0 array_size 0 array_size 0 array_size])
    title(['Iter = ',num2str(iter)])
    grid on
    pause(0.1)
    drawnow
    
    X = mod(X + deltaX, array_size);
    Y = mod(Y + deltaY, array_size);
    Z = mod(Z + deltaZ, array_size);
end



