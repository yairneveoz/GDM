% Charge_Distribution_Atmosphere

clear
clc

data = readmatrix('Charge_Distribution_Atmosphere.csv');

z = data(:,1);
charge = data(:,2);

figure(11)
subplot(1,2,1)
plot(z, charge, '.-')
xlabel('Z (km)')
ylabel('Charge')

subplot(1,2,2)
plot(charge, z, '.-')
xlabel('Charge')
ylabel('Z (km)')

V_0 = (charge(2) - charge(1))/(z(2) - z(1));



