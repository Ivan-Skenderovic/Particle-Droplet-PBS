clear;
clc;
format short e;
%close all;
fclose all;

% define mesh
noSpatialGridpoints = 50000;
x = linspace(0, 1, noSpatialGridpoints);
t = linspace(0, 1e-1, 100);
options = odeset('relTol', 1e-6, 'absTol', 1e-6);

m = 2; % m=2 for spherical coordinates
cb = pdepe(m,@heat,@heatic,@heatbc,x,t,options); % run solver
%plot across a secton for different times
% for i = 1:length(t)
%   col=[i/length(t) 0 1-i/length(t)];
%   plot(x,cb(i,:),'color',col,'linewidth',2); hold on; grid on; box on;
% end
% xlabel('r/r_0 / -')
% ylabel('Temperature / K')

%calculate using method of lines
initialTemperatures(1:noSpatialGridpoints) = 300;
odesystemhandle = @(t,Y) HeatDiffusion_1D_convectionBC(t, Y);
[times, Temperatures] = ode15s(odesystemhandle, t, initialTemperatures, options);

%compare results from pdepe and method of lines
temps_pdepe = cb(end, :);
temps_mol = Temperatures(end, :);
dev = abs((temps_pdepe - temps_mol)./temps_pdepe);

figure(3)
hold on; 
grid on; 
box on;
plot(x, dev);

%%FUNCTION DEFINITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = heat(x,t,u,dudx) % heat diffusion equation
D = 1; % diffusion coefficient
c = 1;
f = D*dudx;
s = 0;
end

function u0 = heatic(x) % initial condition
u0 = 300; %in Kelvin, initially at standard temperature
end

function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t) %BCs
D = 1; % diffusion coefficient
pl = 0;
ql = -D; %symmetry condition at center
h = 1; %heatflux
pr = (ur - 400)*h; %  h*(T - T_inf), convective heating at surface
qr = D; %
end
%%FUNCTION DEFINITIONS END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
