% Space_density_simulation
clear
clc

array_size = 200;
N = 1000;

X0 = array_size*rand(N,1);
Y0 = array_size*rand(N,1);

%
%% 'Lenard Jones'
rr2 = 0.5:0.1:100;

V0 = 2; %10, 1;
scale_r = 10;
shift_up = 0.001; % 0.03
% VLJ = V0*((scale_r./rr2).^4 - (scale_r./rr2).^2 + shift_up);
% VLJ = VR42(V0,scale_r,shift_up,rr2);
% FLJ = FR53(V0,scale_r,rr2) + 0.02;
VR = VR21(V0,scale_r,shift_up,rr2);
% FR = FR32(V0,scale_r,rr2) + 0.0;

f1 = 5;
f2 = -2;
f3 = 0.001;
r1 = 3;
r2 = 10;
% FR = FRsteps(f1,f2,f3,r1,r2,rr2);
FR = FRr2(V0,r1,scale_r,shift_up,rr2);
figure(5)
% plot(rr2, VR, '-')

plot(rr2, FR, '-','LineWidth',2)
hold on
plot([rr2(1), rr2(end)], [0 0], 'k--')
hold off
axis([0 20 -10 10])
% legend('V','F','0')
%
%% plot
grad_scale = 0.1;
X = X0;
Y = Y0;
N_iterations = 50;
steps1 = [5,-2,1e-3,3,10];
steps2 = [2,-2,1e-3,3,10];
steps3 = [5,-2,1e-2,3,10];
steps4 = [5,-2,1e-4,3,10];
steps5 = [5,-2,1e-1,3,10];
steps6 = [5,-1,1e-3,3,10];
steps7 = [5,-1,1e-3,3,15];
steps8 = [5,-1,0,3,10];
steps9 = [5,-1,0,3,15];
s = steps8;

figure(1)
subplot(1,3,1)
FR = FRr2(V0,r1,scale_r,shift_up,rr2);
% FR = FRsteps(s(1),s(2),s(3),s(4),s(5),rr2);
plot(rr2, FR, '-','LineWidth',2)
hold on
plot([rr2(1), rr2(end)], [0 0], 'k--')
hold off
axis square
axis([0 50 -10 10])
xlabel('r')
ylabel('F(r)')
title('F(r)')
text(1,s(1)+1,num2str(s(1)))
text(s(4)+1,s(2)-1,num2str(s(2)))
text(s(5)+1,s(3)+1,num2str(s(3)))

for iter = 1:N_iterations
    dist_X = distD(X);
    dist_Y = distD(Y);
    dist_R = sqrt(dist_X.^2 + dist_Y.^2);
    if iter == 1
        subplot(1,3,2)
        plot(X, Y, '.')
        axis equal
        axis([0 array_size 0 array_size])
        title(['Iter = ',num2str(iter)])
        drawnow
    end

%     F_array = FRsteps(s(1),s(2),s(3),s(4),s(5),dist_R);
    F_array = FRr2(V0,r1,scale_r,shift_up,dist_R);

    deltaX_array = F_array.*dist_X./dist_R;
    deltaY_array = F_array.*dist_Y./dist_R;
    deltaX_array(isnan(deltaX_array)) = 0;
    deltaY_array(isnan(deltaY_array)) = 0;
    deltaX = sum(deltaX_array,2)*grad_scale;
    deltaY = sum(deltaY_array,2)*grad_scale;

    subplot(1,3,3)
    plot(X, Y, '.')
%     hold on
%     quiver(X, Y, deltaX, deltaY) 
%     hold off
    axis equal
    axis([0 array_size 0 array_size])
    title(['Iter = ',num2str(iter)])
    pause(0.1)
    drawnow
    
    X = mod(X + deltaX, array_size);
    Y = mod(Y + deltaY, array_size);
end



