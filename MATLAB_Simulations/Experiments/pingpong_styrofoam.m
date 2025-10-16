% pingpong_styrofoam
clear
clc


roh_styrofoam1 = 0.03; % gr/cm^3. from chatGPT.
m_styrofoam_measured = 1.73; % gr/cm^3. from my measurement.

r_pingpong = 2; % cm
r_styrofoam1 = 2.8; % cm
r_styrofoam3 = 1.2; % cm
r_styrofoam4 = 0.1; % cm

m_pingpong = 2.8; % gr

v_pingpong = volumeR(r_pingpong);
s_pingpong = surfaceR(r_pingpong);
sm_ratio_pingpong = s_pingpong/m_pingpong;

% 
v_styrofoam1 = volumeR(r_styrofoam1);
v_styrofoam3 = volumeR(r_styrofoam3);
v_styrofoam4 = volumeR(r_styrofoam4);

% roh from M and R:
rho_styrofoam2 = m_styrofoam_measured/v_styrofoam1;
rho_styrofoam = rho_styrofoam2;
%
s_styrofoam1 = surfaceR(r_styrofoam1);
s_styrofoam3 = surfaceR(r_styrofoam3);
s_styrofoam4 = surfaceR(r_styrofoam4);

m1_styrofoam = v_styrofoam1*rho_styrofoam;
m2_styrofoam = (v_styrofoam1 - v_pingpong)*rho_styrofoam;
m3_styrofoam = v_styrofoam3*rho_styrofoam;
m4_styrofoam = v_styrofoam4*rho_styrofoam;

sm_ratio_styrofoam1 = s_styrofoam1/m1_styrofoam;
sm_ratio_styrofoam2 = s_styrofoam1/m2_styrofoam;
sm_ratio_styrofoam3 = s_styrofoam3/m3_styrofoam;
sm_ratio_styrofoam4 = s_styrofoam4/m4_styrofoam;

disp(['S to m Pingpong = ', num2str(sm_ratio_pingpong)])
disp(['S to m Full 5.6 Styrofoam = ',num2str(sm_ratio_styrofoam1)])
disp(['S to m Hollow 5.6 Styrofoam = ',num2str(sm_ratio_styrofoam2)])
disp(['S to m Full 2.4 Styrofoam = ',num2str(sm_ratio_styrofoam3)])
disp(['S to m Full 0.2 Styrofoam = ',num2str(sm_ratio_styrofoam4)])

function surface_r = surfaceR(R)
    surface_r = 4*pi*R^2;
end

function volume_r = volumeR(R)
    volume_r = (4/3)*pi*R^3;
end




