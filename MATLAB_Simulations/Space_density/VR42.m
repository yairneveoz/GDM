function V_r = VR42(V0,scale_r,shift_up,dist_R)

V_r = V0*((scale_r./dist_R).^4 - (scale_r./dist_R).^2 + shift_up);

end