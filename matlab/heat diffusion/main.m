%calculate droplet temperature profile during combustion of a pure ethanol
%droplet as system of ordinary differential equations. The surface 
%temperature is then used to calculate droplet evaporations rate
%using ranz and marshall correction
clear;
clc;
format long e;

%materialproperties
NOSPATIALGRIDPOINTS = 100;
DROPLET_INITIAL_TEMPERATURE = 300;

LiqData = LiquidPhaseProperties('Octane');
%LiqData = LiquidPhaseProperties('Ethanol');

GasData = GasPhaseProperties;

%benchmark with C.K. Law paper 

r_0 = 25e-6;
delta_r = r_0/NOSPATIALGRIDPOINTS;

eventFcnHandle = @(t, Y) dropletEvents(t, Y);

options = odeset('relTol', 1e-9, 'absTol', 1e-9, 'Events', eventFcnHandle);

tspan = logspace(-10, 1, 1000);

droplet_initialTemperaturesAndRsquare(1:NOSPATIALGRIDPOINTS) = DROPLET_INITIAL_TEMPERATURE;
droplet_initialTemperaturesAndRsquare(NOSPATIALGRIDPOINTS + 1) = r_0^2;

odesystemhandle = @(t,Y) HeatDiffusion_1D(t, Y, LiqData, GasData);
%odesystemhandle = @(t,Y) HeatDiffusion_1D_uniform_r(t, Y, boundaryType, LiqData, GasData);
tic
[times, TandRsquare] = ode15s(odesystemhandle, tspan, droplet_initialTemperaturesAndRsquare, options);
toc

%%
r = linspace(0, 1, NOSPATIALGRIDPOINTS);

figure(1)
hold on
box on
plot(r, TandRsquare(1, 1:end-1 ),'Color','k'); %matlab 
plot(r, TandRsquare(end, 1:end-1 ),'Color','k'); %matlab 
xlabel('r/r_0');
ylabel('Temperature [K]');
xlim([0 1]);
r_matlab = sqrt(TandRsquare(:,end));

figure(2);
hold on;
box on;
plot(times, r_matlab.^2./r_0.^2,'Color','k');
xlabel('t');
ylabel('R^2');

surfaceTemperature = TandRsquare(:,NOSPATIALGRIDPOINTS);
calcTransportProperties

figure(3);
hold on;
box on;
plot(times, m_d,'Color','k');
xlabel('t');
ylabel('m_d');

figure(4);
hold on;
box on;
plot(times, H_d,'Color','k');
xlabel('t');
ylabel('H_d');




