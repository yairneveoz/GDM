% test_voronoi3

clear
clc

blue_white_red = blueWhiteRedColormap(64);

x = [2 5; 2 5; 8 8];
y = [4 0; 8 2; 4 0];
c = [0; 1];
figure(6)
patch(x,y,c)
colormap(blue_white_red)
colorbar

%%
x = gallery('uniformdata',[10 2],5);
[v,c] = voronoin(x);
for i = 1:length(c)
  if all(c{i}~=1)
    x = v(c{i},1);
    y = v(c{i},2);
    a = polyarea(x,y);
    patch(x,y,a);
  end
end
%