clear;
clc;
format short e;
%close all;
fclose all;

%space and time discretization
NOSPATIALGRIDPOINTS = 1000;
r = linspace(0, 1, NOSPATIALGRIDPOINTS);
TIME_START = 0;
TIME_END = 0.1;
NOTIMEPOINTS = 100;
time = linspace(TIME_START, TIME_END, NOTIMEPOINTS);
PDEPE_ERROR_TOLERANCE = 1e-6;
options = odeset('relTol', PDEPE_ERROR_TOLERANCE, 'absTol', PDEPE_ERROR_TOLERANCE);

m = 2; % m=2 for spherical coordinates
results = pdepe(m, @HeatDiffusion, @HeatDiffusionIC, @HeatDiffusionBC, ...
    r, time,options); % run solver

%calculate using method of lines
initialTemperatures(1:NOSPATIALGRIDPOINTS) = 300;
odesystemhandle = @(t,Y) HeatDiffusion_MOL(t, Y);
[times, Temperatures] = ode15s(odesystemhandle, time, initialTemperatures, options);

%compare results from pdepe and method of lines
temperatures_pdepe = results(end, :);
temperatures_MOL = Temperatures(end, :);
dev = abs((temperatures_pdepe - temperatures_MOL)./temperatures_pdepe)*100;

figure(1)
hold on; 
grid on; 
box on;
plot(r, temperatures_pdepe, 'LineStyle','-');
plot(r, temperatures_MOL, 'LineStyle','--');
xlabel('Radial coordinate r/r_0 / -')
ylabel('Temperature / K')

%plot relative deviation between method of lines and pdepe
figure(2)
hold on; 
grid on; 
box on;
plot(r, dev);
xlabel('Radial coordinate r/r_0 / -')
ylabel('Rel. deviation / % ')

function [c,f,s] = HeatDiffusion(r, time, u, dTdr) 
    D = 1; %diffusion coefficient
    c = 1;
    f = D*dTdr;
    s = 0;
end

function u0 = HeatDiffusionIC(r) % initial condition
    u0 = 300; %in Kelvin, initially at standard temperature
end

function [pl,ql,pr,qr] = HeatDiffusionBC(rl, ul, rr, ur, time) 
    %define boundary conditions for pdepe
    D = 1; %diffusion coefficient
    pl = 0;
    ql = -D; %symmetry condition at center
    h = 1; %heatflux
    pr = (ur - 400)*h; %
    qr = D; %
end

