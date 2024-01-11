clear;
clc;
format long e;

NOSPATIALGRIDPOINTS = 100;
DROPLET_INITIAL_TEMPERATURE = 300; %K
DROPLET_INITIAL_RADIUS = 100e-6;

%material properties
LiqData = LiquidPhaseProperties('Octane');
GasData = GasPhaseProperties;

%setup ode solver
ODEINT_ERROR_TOLERANCE = 1e-9;
odesystemhandle = @(t,Y) HeatDiffusion_1D(t, Y, LiqData, GasData);
options = odeset('relTol', ODEINT_ERROR_TOLERANCE, ...
    'absTol', ODEINT_ERROR_TOLERANCE);

TIME_START = 1e-10; %s
TIME_END = 1; %s
NOTIMESTEPS = 1000;

tspan = logspace(log10(TIME_START), log10(TIME_END), NOTIMESTEPS);
droplet_initialTemperaturesAndRsquare(1:NOSPATIALGRIDPOINTS) = DROPLET_INITIAL_TEMPERATURE;
droplet_initialTemperaturesAndRsquare(NOSPATIALGRIDPOINTS + 1) = DROPLET_INITIAL_RADIUS^2;

[times, TandRsquare] = ode15s(odesystemhandle, tspan, droplet_initialTemperaturesAndRsquare, options);

%normalized radial coordinate
r = linspace(0, 1, NOSPATIALGRIDPOINTS);
Rsquare = TandRsquare(:,end);

%plot temperatures at start and end of simulation
figure(1)
hold on
box on
plot(r, TandRsquare(1, 1:end-1 ),'Color','k');
plot(r, TandRsquare(end, 1:end-1 ),'Color','k');
xlabel('r/r_0');
ylabel('Temperature [K]');
xlim([0 1]);

%plot normalized droplet surface over time
figure(2);
hold on;
box on;
plot(times, Rsquare./DROPLET_INITIAL_RADIUS^2, 'Color','k');
xlabel('t');
ylabel('R^2/R^2_0');

