% test_voronoi4
clear
clc

blue_white_red = flipud(blueWhiteRedColormap(64));

rng default
figure(7)
pts = randn(500,2);
[v,c] = voronoin(pts);
for i = 1:length(c)
  if all(c{i}~=1)
    x = v(c{i},1);
    y = v(c{i},2);
    a = polyarea(x,y);
    patch(x,y,a);
  end
end
hold on
scatter(pts(:,1),pts(:,2),40,[.9 .1 .1],'Marker','.','MarkerEdgeAlpha',.5)
colormap(blue_white_red)
colorbar('SouthOutside')
axis equal
xlim([-3 3])
ylim([-3 3])
caxis([0 1])


