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

CD = CD0 + (1 - sigma) .* K .* (CL0 + CLalpha.*x).^2 + sigma .* (2.*sign(x).*sin(x).^3);
CDreg = CD0 + K .* (CL0 + CLalpha .* x) .^2;
CDdsl = 2*sign(x).*(sin(x).^3);

CM = CM0 + CMalpha * x;

% L = 1/2 * density_air .* V.^2 .* (CL + ((CLq .* mean_chord.*q) ./ (2.*V)) .* CLdelta.*delta);
% D = 1/2 * density_air .* V.^2 .* (CD + ((CDq .* mean_chord.*q) ./ (2.*V)) .* CDdelta.*delta);
% M = 1/2 * density_air .* V.^2 .* mean_chord .* (CM + ((CMq .* mean_chord.*q) ./ (2.*V)) .* CMdelta.*delta);