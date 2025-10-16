function F_r = FR53(V0,scale_r,dist_R)

F_r = -V0*(4*(scale_r^4./dist_R.^5) - 2*(scale_r^2./dist_R.^3));

end