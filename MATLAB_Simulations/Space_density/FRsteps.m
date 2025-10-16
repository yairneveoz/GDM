function F_r = FRsteps(f1,f2,f3,r1,r2,dist_R)

F_r = zeros(size(dist_R));
F_r(dist_R < r1) = f1;
F_r(dist_R >= r1) = f2;
F_r(dist_R >= r2) = f3;

end