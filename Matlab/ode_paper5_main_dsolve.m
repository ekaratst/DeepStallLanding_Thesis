syms v(t) y(t) q(t) o(t) h(t) r(t)
Iyy = 0.388;
m = 1.5;
g = 9.81;
density_air = 1.225;
wing_area = 0.28;
mean_chord = 0.31;
CL0 = 0.062;
CLalpha = 6.098;
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
T = 0;
%<<<-----------------angle----------------->>>
theta = 30;
gramma = -20;
aoa = theta - gramma;
% delta = [0,-25,-50,-75];
% delta = [0,-15,-45,-70];  %useeee
elevator_angle = -5;
pitch_angle = theta*pi/180;
angle_of_descent = gramma*pi/180; 
x = aoa*pi/180;
delta = elevator_angle*pi/180;
%<<<--------------------------------------->>>
aircraft_velocity = 8;
pitch_rate = 0;
height = 0;
horizontal_distance = 0;
%-----------------------equations-------------------------------------------------------------------
sigma = ((1 + exp(-M*(x - alpha0)) + exp(M*(x + alpha0))) ./ ((1 + exp(-M*(x - alpha0))) .* (1 + exp(M*(x + alpha0)))));
CL = (1 - sigma).*(CL0 + CLalpha.*x) + sigma.*(2.*sign(x).*sin(x).^2 .* cos(x));
CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
CM = CM0 + CMalpha * x;
t_interval = [0,10];

ode1 = diff(v) == ((T*cos(x)) - (1/2 * density_air * wing_area * v^2 * (CD + ((CDq * mean_chord*q) / (2*v)) + CDdelta*delta)) - (m*g*sin(y))) / m; 
ode2 = diff(y) == ((T*sin(x)) + (1/2 * density_air * wing_area * v^2 * (CL + ((CLq * mean_chord*q) / (2*v)) + CLdelta*delta)) - (m*g*cos(y))) / m*v; 
ode3 = diff(q) == (1/2 * density_air * wing_area * v^2 * mean_chord * (CM + ((CMq * mean_chord*q) / (2*v)) + CMdelta*delta)) / Iyy; 
ode4 = diff(o) == q; 
ode5 = diff(h) == v*sin(y); 
ode6 = diff(r) == v*cos(y);

odes = [ode1; ode2; ode3; ode4; ode5; ode6];
S  = dsolve(odes);
vSol(t) = S.v;
ySol(t) = S.y;
qSol(t) = S.q;
oSol(t) = S.o;
hSol(t) = S.h;
rSol(t) = S.r;
disp(vSol(t));

