function V_r = VR21(V0,scale_r,shift_up,dist_R)

V_r = V0*((scale_r./dist_R).^2 - (scale_r./dist_R).^1 + shift_up);

end