% Spirals2
clear
clc

R0 = 137.0366;

loops = 1;
c = 1;
v1 = 0.0*c;
v2 = 0.1*c;
v3 = 0.5*c;
v4 = 0.8*c;
v5 = 0.95*c;

c_alpha = 0.5;
N_samples = 5;
cmap = jet(N_samples);

[R_path1, r_path1] = getSpiralPath(v1,loops);
[R_path2, r_path2] = getSpiralPath(v2,loops);
[R_path3, r_path3] = getSpiralPath(v3,loops);
[R_path4, r_path4] = getSpiralPath(v4,loops);
[R_path5, r_path5] = getSpiralPath(v5,loops);

figure(14)
clf
plot3(r_path1(:,1), r_path1(:,2), r_path1(:,3), 'b-','LineWidth',1)
hold on
plot3(R_path1(:,1), R_path1(:,2), R_path1(:,3), '-','LineWidth',3,...
    'Color',[cmap(1,:), c_alpha])

plot3(r_path2(:,1), r_path2(:,2), r_path2(:,3), 'b-','LineWidth',1)
plot3(R_path2(:,1), R_path2(:,2), R_path2(:,3), '-','LineWidth',3,...
    'Color',[cmap(2,:), c_alpha])

% plot3(r_path3(:,1), r_path3(:,2), r_path3(:,3), 'b-','LineWidth',1)
% plot3(R_path3(:,1), R_path3(:,2), R_path3(:,3), '-','LineWidth',3,...
%     'Color',[cmap(3,:), c_alpha])

% plot3(r_path4(:,1), r_path4(:,2), r_path4(:,3), 'b-','LineWidth',1)
% plot3(R_path4(:,1), R_path4(:,2), R_path4(:,3), '-','LineWidth',3,...
%     'Color',[cmap(4,:), c_alpha])

% plot3(r_path5(:,1), r_path5(:,2), r_path5(:,3), 'b-','LineWidth',1)
% plot3(R_path5(:,1), R_path5(:,2), R_path5(:,3), '-','LineWidth',3,...
%     'Color',[cmap(5,:), c_alpha])
hold off

xlabel('X(r_0)')
ylabel('Y(r_0)')
zlabel('Z(r_0)')

axis equal
% axis(1.2*[-R0 R0 -R0 R0 -1 max(Z)])
grid on
