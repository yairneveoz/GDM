% test_voronoi2
X = [-1.5 3.2; 1.8 3.3; -3.7 1.5; -1.5 1.3; ...
0.8 1.2; 3.3 1.5; -4.0 -1.0;-2.3 -0.7; ...
0 -0.5; 2.0 -1.5; 3.7 -0.8; -3.5 -2.9; ...
-0.9 -3.9; 2.0 -3.5; 3.5 -2.25];

X(:,3) = [ 1 2 1 3 1 2 2 2 2 3 3 3 3 3 3]';
%% %
blue_white_red = blueWhiteRedColormap(64);

% figure(5)
% pcolor(rand(20))
% colormap(blue_white_red)
%% %
ccode = ["red","green","blue"];

dt = delaunayTriangulation(X(:,1:2));
[V,R] = voronoiDiagram(dt);
figure(4)
voronoi(X(:,1),X(:,2))
hold on
for i = 1:size(X,1)
    A=V(R{i},:);
    B=A(any(~isinf(A),2),:);
    if(size(B,1)>2)
        plot(polyshape(B),'FaceColor',ccode(X(i,3)));
    end
end


