function F_r = FRr2(V0,r1,scale_r,shift_up,dist_R)

F_r = -V0*(scale_r./dist_R.^2 + shift_up);
F_r(dist_R < r1) = V0;

end