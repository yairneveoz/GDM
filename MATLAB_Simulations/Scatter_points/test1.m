close all; 
clc;

N = 10;
x = rand(N,1);
y = rand(N,1);

% 1) Make sure we're calling the built-in functions
% disp(which('plot','-all'));
% disp(which('scatter','-all'));

figure('Color','w','Position',[100 100 640 480],'Visible','on');
ax = axes('Parent', gcf); 
grid(ax,'on'); 
box(ax,'on'); 
hold(ax,'off');

% 2) Try a very explicit low-level marker (bypasses some plot defaults)
h = line(ax, x, y, 'LineStyle','none', 'Marker','.', 'MarkerSize',20, 'Color','b');
title(ax, 'Low-level line() with markers');
xlabel(ax,'X'); ylabel(ax,'Y');
drawnow;
fprintf('Children in axes after line(): %d\n', numel(ax.Children));

% pause(1);
% 
% % 3) If still blank, switch renderer and try scatter
% set(gcf,'Renderer','painters');  % safest 2D renderer
% clf; ax = axes('Parent', gcf); grid(ax,'on'); box(ax,'on');
% scatter(ax, x, y, 80, 'r', 'filled');
% title(ax, 'scatter() using painters');
% xlabel(ax,'X'); ylabel(ax,'Y');
% drawnow;
% fprintf('Renderer now: %s\n', get(gcf,'Renderer'));
% fprintf('Children in axes after scatter(): %d\n', numel(ax.Children));

% pause(1);
% 
% % 4) As a last resort, force software OpenGL and try plot
% opengl software;
% clf; ax = axes('Parent', gcf); grid(ax,'on'); box(ax,'on');
% plot(ax, x, y, 'b.', 'MarkerSize', 20);
% title(ax, 'plot() after opengl software');
% xlabel(ax,'X'); ylabel(ax,'Y');
% drawnow;
% fprintf('OpenGL mode forced to software.\n');
