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
delta = 20*pi/180;

% x_degrees = linspace(-10,180);
x_degrees = 30;
x = x_degrees.*pi./180;
%x_degrees = [-10 -5 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];

sigma = ((1 + exp(-M*(x - alpha0)) + exp(M*(x + alpha0))) ./ ((1 + exp(-M*(x - alpha0))) .* (1 + exp(M*(x + alpha0)))));

CL = (1 - sigma).*(CL0 + CLalpha.*x) + sigma.*(2.*sign(x).*sin(x).^2 .* cos(x));
CLreg = CL0 + CLalpha .* x;
CLdsl = 2*sign(x).*(sin(x).^2) .* cos(x);

CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
CDreg = CD0 + K .* (CL0 + CLalpha .* x) .^2;
CDdsl = 2*sign(x).*(sin(x).^3);

CM = CM0 + CMalpha * x;

% L = 1/2 * density_air * S .* V.^2 .* (CL + ((CLq .* mean_chord.*q) ./ (2.*V)) .* CLdelta.*delta);
% D = 1/2 * density_air * S .* V.^2 .* (CD + ((CDq .* mean_chord.*q) ./ (2.*V)) .* CDdelta.*delta);
% M = 1/2 * density_air * S .* V.^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*q) ./ (2.*V)) .* CMdelta.*delta);

% syms t x a
% g = @(t,x,a)[-x(1)+a*x(3);-x(2)+2*x(3);x(1)^2-2*x(3)];
% for a = 0:2
% [t,xa] = ode45(@(t,x) g(t,x,a),[0 1.5],[1 1/2 3]);
% figure
% plot(t,xa(:,2))
% title(['y(t) for a=',num2str(a)'])
% end


% syms t x
% L = (1/2 * density_air * S .* x(1).^2 .* (CL + ((CLq .* mean_chord.*x(3)) ./ (2.*x(1))) .* CLdelta.*delta));
% D = (1/2 * density_air * S .* x(1).^2 .* (CD + ((CDq .* mean_chord.*x(3)) ./ (2.*x(1))) .* CDdelta.*delta));
% M = (1/2 * density_air * S .* x(1).^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*x(3)) ./ (2.*x(1))) .* CMdelta.*delta));
% equation_of_motion = @(t,x)[(-(1/2 * density_air * S .* x(1).^2 .* (CD + ((CDq .* mean_chord.*x(3)) ./ (2.*x(1))) .* CDdelta.*delta))) - (m*g*sin(x(2))); 
%     (1/2 * density_air * S .* x(1).^2 .* (CL + ((CLq .* mean_chord.*x(3)) ./ (2.*x(1))) .* CLdelta.*delta)) - (m*g*sin(x(2))); 
%     (1/2 * density_air * S .* x(1).^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*x(3)) ./ (2.*x(1))) .* CMdelta.*delta))./Iyy; 
%     x(3); 
%     x(1)*sin(x(2)); 
%     x(1)*cos(x(2))];
% equation_of_motion = @(t,x)[(-D) - (m*g*sin(x(2))); 
%     L - (m*g*sin(x(2))); 
%     M./Iyy; 
%     x(3); 
%     x(1)*sin(x(2)); 
%     x(1)*cos(x(2))];
% [t,x1] = ode45(@(t,x) equation_of_motion(t,x), [0 2], [1,2,3]);
% figure
% plot(t,x1(:,2))
% title(['y(t)'])



%time interval and initial conditions
t_interval = [0,1000];
init_cond = [10,5*pi/180,0,10*pi/180,15,30]';
%solution
[t,y] = ode45(@(t,Y) odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy) , t_interval , init_cond);
%plot
plot(t,y(:,1),'b',t,y(:,2),'r');


function dYdt = odefcn(t,Y,density_air,S,CL,CLq,mean_chord,CLdelta,delta,CD,CDq,CDdelta,CM,CMq,CMdelta,m,g,Iyy)
dYdt = [    (-(1/2 * density_air * S .* Y(1).^2 .* (CD + ((CDq .* mean_chord.*Y(3)) ./ (2.*Y(1))) .* CDdelta.*delta))) - (m*g*sin(Y(2))); 
            (1/2 * density_air * S .* Y(1).^2 .* (CL + ((CLq .* mean_chord.*Y(3)) ./ (2.*Y(1))) .* CLdelta.*delta)) - (m*g*sin(Y(2))); 
            (1/2 * density_air * S .* Y(1).^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*Y(3)) ./ (2.*Y(1))) .* CMdelta.*delta))./Iyy; 
            Y(3); 
            Y(1)*sin(Y(2)); 
            Y(1)*cos(Y(2))];
end

