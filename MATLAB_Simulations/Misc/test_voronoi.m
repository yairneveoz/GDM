% test_voronoi
clear
clc

X = rand(20,1);
Y = rand(20,1);
P = [X,Y];
% [V, C] = voronoin(P);

[v, c] = voronoin(P);
cmap = [
    1 0 0;  % red
    0 1 0;  % green
    0 0 1;  % blue
    1 1 0];  % yellow

h = patch('XData', v(:,1), 'YData', v(:,2), 'CData', cmap(c,:));



