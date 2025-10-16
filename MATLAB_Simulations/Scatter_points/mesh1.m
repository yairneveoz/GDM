% Define the number of points along one dimension

clear
clc
close all

% pdetool
% pdeplot(p,e,t,'XYData',u,'Contour','on'); 
% axis equal tight
% % Electric field magnitude (finite-difference approx):
% [ux,uy] = pdegrad(p,t,u);
% Emag = sqrt(ux.^2 + uy.^2);  % at triangle centers

model = createpde();

% Geometry: rectangle 0<=x<=0.2, 0<=y<=0.1
R1 = [3 4 0 0.2 0.2 0 ...
          0 0   0.1 0.1]';
gd = R1; ns = char('R1')'; sf = 'R1';
g = decsg(gd,sf,ns);
geometryFromEdges(model,g);

% Coefficients: Laplace (c=1, a=0, f=0), steady (d=0)
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0);

% Boundaries: left x=0 -> u=+1 ; right x=0.2 -> u=-1
applyBoundaryCondition(model,'dirichlet','Edge',1,'u',1);   % edge IDs may differ; use pdegplot(model,'EdgeLabels','on')
applyBoundaryCondition(model,'dirichlet','Edge',3,'u',-1);

% Mesh & solve
generateMesh(model,'Hmax',0.005);
result = solvepde(model);
u = result.NodalSolution;

pdeplot(model,'XYData',u,'Contour','on'); axis equal tight
