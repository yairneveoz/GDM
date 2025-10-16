% random_points_dynamics_1D

clear
clc

N = 200; % number of points
L = 10; % smaple length

x0 = sort(L*rand(N,1));
x = x0;

L_factor = 1;
T_factor = 1;
dx0 = L/N;
% DX0 = dx0*ones(size(x));
DX0 = dx0*0.5*randn(size(x));
figure(3)
subplot(3,1,1)
plot(x0, 1*ones(size(x0)), '.')
ylim([0 2])
xlabel('x')
title('Initial spatial distribution')

for iter = 1:1000
    disp(iter)
    %%%
    mean_diff_x = zeros(size(x));
    for n = 2:N-1
%         mean_diff_x(n) = mean([abs(x(n) - x(n-1)),...
%             abs(x(n) - x(n+1))]);
        mean_diff_x(n) = mean([abs(x(n) - x(n-1)),...
            abs(x(n) - x(n+1))]);
    end
    mean_diff_x(1) = abs(x(2) - x(1));
    mean_diff_x(N) = abs(x(N) - x(N-1));
    %%%

    subplot(3,1,2)
    plot(mean_diff_x, '.-')
%     plot(mean_diff_x/dx0, '.-')
    % ylim([0 2])
    xlabel('Cell index')
    ylabel('Cell size')
    title('Cell sizes by index')
    
    DX = DX0.*mean_diff_x;
%     x1 = x + DX;
    %%% Periodic
    x1_0 = mod(x + DX, L);
    x1 = sort(x1_0);
    %%%

    subplot(3,1,3)
    plot(x1, 1*ones(size(x0)), 'r.')
    ylim([0 2])
    pause(0.05)
    drawnow
    
    x = x1;
    
end