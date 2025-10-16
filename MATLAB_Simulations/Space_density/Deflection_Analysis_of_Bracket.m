% Deflection_Analysis_of_Bracket

clear
clc

model = createpde('structural','static-solid');
importGeometry(model,'BracketWithHole.stl');

figure(11)
clf
subplot(2,4,1)
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Bracket with Face Labels')
%
subplot(2,4,2)
pdegplot(model,'FaceLabels','on')
view(-134,-32)
title('Bracket with Face Labels, Rear View')

structuralProperties(model,'YoungsModulus',200e9, ...
                           'PoissonsRatio',0.3);
structuralBC(model,'Face',4,'Constraint','fixed');
structuralBoundaryLoad (model,'Face',8,'SurfaceTraction',...
    [0;0;-1e4]);
generateMesh(model);
%
subplot(2,4,3)
pdeplot3D(model)
title('Mesh with Quadratic Tetrahedral Elements');

result = solve(model);
minUz = min(result.Displacement.uz);
fprintf('Maximal deflection in the z-direction is %g meters.',...
    minUz)

subplot(2,4,5)
pdeplot3D(model,'ColorMapData',result.Displacement.ux)
title('x-displacement')
colormap('jet')


subplot(2,4,6)
pdeplot3D(model,'ColorMapData',result.Displacement.uy)
title('y-displacement')
colormap('jet')

subplot(2,4,7)
pdeplot3D(model,'ColorMapData',result.Displacement.uz)
title('z-displacement')
colormap('jet')

subplot(2,4,8)
pdeplot3D(model,'ColorMapData',result.VonMisesStress)
title('von Mises stress')
colormap('jet')



