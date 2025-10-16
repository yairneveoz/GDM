function F_r = FR32(V0,scale_r,dist_R)

f2_scale = 2; % 1
F_r = -V0*(2*(scale_r^2./dist_R.^3) - f2_scale*(scale_r./dist_R.^2));
F_r(F_r > V0) = V0;

end