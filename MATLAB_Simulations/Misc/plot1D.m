function [] = plot1D(inputs,outputs)


rho_0 = inputs.rho_0;

r0 = outputs.r0;
rho1_p = outputs.rho1_p;
rho1_n = outputs.rho1_n;
outputs.rho1_p = rho1_p;
outputs.rho1_n = rho1_n;
rp = outputs.rp;
rn = outputs.rn;
q1p = outputs.q1p;
q1n = outputs.q1n;
Q1p = outputs.Q1p;
Q1n = outputs.Q1n;
E1p = outputs.E1p;
E1n = outputs.E1n;
m1p = outputs.m1p;
m1n = outputs.m1n;
half_sum_q10p = outputs.half_sum_q10p;
half_sum_q10n = outputs.half_sum_q10n;

figure(19)
clf

%% (2,4,1):
subplot(2,4,1)
plot(r0,rho1_p,'r.-')
hold on
plot(r0,rho_0*ones(size(rho1_p)),'.-', 'Color', 0.5*[1 1 1])
plot(r0,rho1_n,'b.-')
hold off

xlabel('r_0')
ylabel('\rho(r_0)/\rho_0(r_0)')
title('\rho(r_0)')
legend({'\rho_+','\rho_0','\rho_-'})
% xticks(0:10:floor(grid_size/2))
%
%% (2,4,2):
subplot(2,4,2)
plot(r0, q1p, 'r.-')
hold on
plot(r0, 0*ones(size(r0)), '.-', 'Color',0.5*[1 1 1])
plot(r0, q1n, 'b.-')
plot(4.1*[1 1],[0 max(q1p)],'r--')
plot(4.1*[1 1],[min(q1n) 0],'b--')
hold off
% axis([0 30 0 inf])
xlabel('r_0')
ylabel('q')
title('q(r_0)')
legend({'q_+','q_0','q_-'})
%
%% (2,4,3):
subplot(2,4,3)
plot(r0, Q1p, 'r.-')
hold on
plot(r0, Q1n, 'b.-')
hold off
xlabel('r_0')
ylabel('Q')
title('Q(r_0)')
legend({'Q_+','Q_-'})
%
%% (2,4,4):
subplot(2,4,4)
plot(r0, E1p, 'ro-')
hold on
plot(r0, E1n, 'b.-')
hold off
xlabel('r_0')
ylabel('U')
title('U(r_0)')
legend({'U_+','U_-'})
%
%% (2,4,5):
subplot(2,4,5)
plot(rp, 1.5*ones(size(rp)),'r.-')
hold on
plot(r0, 1.0*ones(size(r0)),'.-', 'Color', 0.5*[1 1 1])
plot(rn, 0.5*ones(size(rn)),'b.-')
hold off
axis([0 30 0 2])
xlabel('r')
% ylabel('\rho/\rho_0')
title('r')
legend({'r_+','r_0','r_-'})
%
%% (2,4,6):
subplot(2,4,6)
plot(rp, q1p, 'r.-')
hold on
plot(r0, 0*ones(size(r0)), '.-', 'Color', 0.5*[1 1 1])
plot(rn, q1n, 'b.-')
plot(0.7*[1 1],[0 max(q1p)],'r--')
plot(7.2*[1 1],[min(q1n) 0],'b--')
hold off
% axis([0 30 0 inf])
xlabel('r')
ylabel('q')
title('q(r)')
legend({'q_+','q_0','q_-'})
%
%% (2,4,7):
subplot(2,4,7)
plot(rp, Q1p, 'r.-')
hold on
plot(rn, Q1n, 'b.-')
hold off
xlabel('r')
ylabel('Q')
title('Q(r)')
legend({'Q_+','Q_-'})
%
%% (2,4,8):
subplot(2,4,8)
plot(rp, E1p, 'r.-')
hold on
plot(rn, E1n, 'b.-')
hold off
xlabel('r')
ylabel('U')
title('U(r)')
legend({'U_+','U_-'})
%
end


