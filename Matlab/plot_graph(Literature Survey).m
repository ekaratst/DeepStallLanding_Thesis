Iyy = 0.388;
m = 1.5;
density_air = 1.225;
S = 0.28;
mean_chord = 1.3;
CL0 = 0.062;
CLalpha = 6.09;
CD0 = 0.098;
K = 0.012;
CM0 = 0.028;
CMalpha = -0.031;
CLq = 0;
CLdelta = -1.72;
CDq = 0;
CDdelta = -0.814;
CMq = -13.1;
CMdelta = -0.325;
alpha0 = 20*pi/180;
M = 50;
delta = 20*pi/180;

x_degrees = linspace(-10,180);
x = x_degrees.*pi./180;
%x_degrees = [-10 -5 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];

sigma = ((1 + exp(-M*(x - alpha0)) + exp(M*(x + alpha0))) ./ ((1 + exp(-M*(x - alpha0))) .* (1 + exp(M*(x + alpha0)))));

CL = (1 - sigma).*(CL0 + CLalpha.*x) + sigma.*(2.*sign(x).*sin(x).^2 .* cos(x));
CLreg = CL0 + CLalpha .* x;
CLdsl = 2*sign(x).*(sin(x).^2) .* cos(x);
sin_test = sin(x).^2 ;

CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
CDreg = CD0 + K .* (CL0 + CLalpha .* x) .^2;
CDdsl = 2*sign(x).*(sin(x).^3);

CM = CM0 + CMalpha * x;

% L = 1/2 * density_air .* V.^2 .* (CL + ((CLq .* mean_chord.*q) ./ (2.*V)) .* CLdelta.*delta);
% D = 1/2 * density_air .* V.^2 .* (CD + ((CDq .* mean_chord.*q) ./ (2.*V)) .* CDdelta.*delta);
% M = 1/2 * density_air .* V.^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*q) ./ (2.*V)) .* CMdelta.*delta);



subplot(2,1,1);
plot(x_degrees, CL);
% plot(x_degrees, sin_test);
hold on
plot(x_degrees, CLreg);
plot(x_degrees, CLdsl);
%xlim([-10 90])
ylim([-2 2])
title('C_L - AOA')
xlabel('Angel of Attack [deg]') 
ylabel('Lift Coefficient') 
legend({'CL','CL_r_e_g','CL_D_S_L'},'Location','northeast')
hold off

subplot(2,1,2);
plot(x_degrees, CD);
hold on
plot(x_degrees, CDreg);
plot(x_degrees, CDdsl)
xlim([-10 90])
%ylim([-5 10])
title('C_D - AOA')
xlabel('Angel of Attack [deg]') 
ylabel('Drag Coefficient') 
legend({'CD','CD_r_e_g','CD_D_S_L'},'Location','northeast')
hold off

% subplot(3,2,2);
% syms V(t) g(t) q(t) theta(t) h(t) r(t)
% ode1 = diff(V) == (-D-m*9.81*sin(g))/m;
% ode2 = diff(g) == (L - m*9.81*cos(g))/(m*V);
% ode3 = diff(q) == M/Iyy;
% ode4 = diff(theta) == q;
% ode5 = diff(h) == V*sin(g);
% ode6 = diff(r) == V*cos(g);
% odes = [ode1; ode2; ode3; ode4; ode5; ode6];
% S = dsolve(odes);
% vSol(t) = S.V;
% gSol(t) = S.g;
% qSol(t) = s.q;
% thetaSol(t) = s.theta;
% hSol(t) = s.h;
% rSol(t) = s.r;
% [vSol(t), gSol(t), qSol(t), thetaSol(t), hSol(t), rSol(t)] = dsolve(odes);
% cond1 = V(0) == 1;
% cond2 = g(0) == 1;
% cond3 = q(0) == 1;
% cond4 = theta(0) == 1;
% cond5 = h(0) == 1;
% cond6 = r(0) == 1;
% conds = [cond1; cond2];
% [vSol(t), gSol(t), qSol(t), thetaSol(t), hSol(t), rSol(t)] = dsolve(odes,conds);
% figure
% title('V - time');
% fplot(vSol)
% grid on
% legend('vSol','Location','best')
% 
