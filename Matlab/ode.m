%parameters
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

% x_degrees = linspace(-10,180);
x_degrees = 30;
x = x_degrees*pi/180;
%x_degrees = [-10 -5 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];

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
init_cond = [8, 5*pi/180, 0, 0, 0, 0]; % [aircraft velocity, angle of descent, pitch angle, pitch rate, height, horizontal distance]

%solution
[t,y] = ode45(@(t,Y) odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy,T,x) , t_interval , init_cond);

% plot(t,y(:,1),'b',t,y(:,2),'r');
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

