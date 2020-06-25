%-----------------------parameters-------------------------------------------------------------------
Iyy = 0.388;
m = 1.5;
g = 9.81;
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
delta = 5*pi/180;
T = 0;

%-------------------------------------AOA----------------------------------------------------------------------------------
x_degrees_aoa = [-10 -5 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];
x_aoa = x_degrees_aoa*pi/180;
sigma_aoa = ((1 + exp(-M*(x_aoa - alpha0)) + exp(M*(x_aoa + alpha0))) ./ ((1 + exp(-M*(x_aoa - alpha0))) .* (1 + exp(M*(x_aoa + alpha0)))));
CL_aoa = (1 - sigma_aoa).*(CL0 + CLalpha.*x_aoa) + sigma_aoa.*(2.*sign(x_aoa).*sin(x_aoa).^2 .* cos(x_aoa));
CLreg_aoa = CL0 + CLalpha .* x_aoa;
CLdsl_aoa = 2*sign(x_aoa).*(sin(x_aoa).^2) .* cos(x_aoa);
CD_aoa = CD0 + (1 - sigma_aoa) .* K .* (CL0 + CLalpha.*x_aoa).^2 + sigma_aoa .* (2.*sign(x_aoa).*sin(x_aoa).^3);
CDreg_aoa = CD0 + K .* (CL0 + CLalpha .* x_aoa) .^2;
CDdsl_aoa = 2*sign(x_aoa).*(sin(x_aoa).^3);
CM_aoa = CM0 + CMalpha * x_aoa;
%-------------------------------------ODE----------------------------------------------------------------------------------
% x_degrees = linspace(-10,180);
x_degrees = 30;
x = x_degrees*pi/180;
sigma = ((1 + exp(-M*(x - alpha0)) + exp(M*(x + alpha0))) ./ ((1 + exp(-M*(x - alpha0))) .* (1 + exp(M*(x + alpha0)))));
CL = (1 - sigma).*(CL0 + CLalpha.*x) + sigma.*(2.*sign(x).*sin(x).^2 .* cos(x));
CLreg = CL0 + CLalpha .* x;
CLdsl = 2*sign(x).*(sin(x).^2) .* cos(x);
CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
CDreg = CD0 + K .* (CL0 + CLalpha .* x) .^2;
CDdsl = 2*sign(x).*(sin(x).^3);
CM = CM0 + CMalpha * x;

%time interval and initial conditions
t_interval = [0,10];

% init_cond = [8, -5*pi/180, -5*pi/180, -15, 0, 0]; 
% init_cond = [8, 0, 0, 0, 0, 0]; 

init_cond = [8,             5*pi/180,          0,            0,        0,           0]; 
% [aircraft velocity | angle of descent | pitch angle | pitch rate | height | horizontal distance]

%solution
[t,y] = ode45(@(t,Y) odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy,T,x) , t_interval , init_cond);
%-------------------------------------AOA-----------------------------------------------------------------------------------
% plot(t,y(:,1),'b',t,y(:,2),'r');
figure(1)
subplot(2,1,1);
plot(x_degrees_aoa, CL_aoa);
% plot(x_degrees_aoa, sin_test);
hold on
plot(x_degrees_aoa, CLreg_aoa);
plot(x_degrees_aoa, CLdsl_aoa);
%xlim([-10 90])
ylim([-2 2])
title('C_L - AOA')
xlabel('Angel of Attack [deg]') 
ylabel('Lift Coefficient') 
legend({'CL','CL_r_e_g','CL_D_S_L'},'Location','northeast')
hold off

subplot(2,1,2);
plot(x_degrees_aoa, CD_aoa);
hold on
plot(x_degrees_aoa, CDreg_aoa);
plot(x_degrees_aoa, CDdsl_aoa)
xlim([-10 90])
%ylim([-5 10])
title('C_D - AOA')
xlabel('Angel of Attack [deg]') 
ylabel('Drag Coefficient') 
legend({'CD','CD_r_e_g','CD_D_S_L'},'Location','northeast')
hold off

%------------------------------------ODE------------------------------------------------------------------------------------

figure(2)
subplot(3,2,1);
plot(y(:,6),y(:,5)); %vetical distance(y-axis)-horizontal distance(x-axis)
xlabel('Horizontal Distance [m]') 
ylabel('Vetical Distance [m]') 
title('Simulated DSL Trajectory')
legend("8 m/s")

subplot(3,2,2);
plot(t,y(:,1))
xlabel('Time') 
ylabel('V [m/s]') 
title('DSL Velocity')
legend("8 m/s")

subplot(3,2,3);
plot(t,y(:,2))
xlabel('Time') 
ylabel('gramma') 
title('DSL Angle of Descent')
legend("8 m/s")

subplot(3,2,4);
plot(t,y(:,3))
xlabel('Time') 
ylabel('theta') 
title('DSL Pitch Angle')
legend("8 m/s")

subplot(3,2,5);
plot(t,y(:,4))
xlabel('Time') 
ylabel('q') 
title('DSL Pitch rate')
legend("8 m/s")

function dYdt = odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy,T,x)
dYdt = [    (T*cos(x)) - (1/2 * density_air * S .* Y(1).^2 .* (CD + ((CDq .* mean_chord.*Y(3)) ./ (2.*Y(1))) + CDdelta.*delta)) - (m*g*sin(Y(2))); 
            (T*sin(x)) + (1/2 * density_air * S .* Y(1).^2 .* (CL + ((CLq .* mean_chord.*Y(3)) ./ (2.*Y(1))) + CLdelta.*delta)) - (m*g*cos(Y(2))); 
            (1/2 * density_air * S .* Y(1).^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*Y(3)) ./ (2.*Y(1))) + CMdelta.*delta)) ./ Iyy; 
            Y(3); 
            Y(1)*sin(Y(2)); 
            Y(1)*cos(Y(2))];
end

