NA = 6.02214086e23; % 1/mol, Avogadro constant

KB = 1.380649e-23; %J/K, Boltzmann constant

NOSTEPS = 1e5;

TIMESTEP = 1e-3; %s

K = 0.8e-6; %droplet surface regression rate in m²/s

ODEINT_ABS_ERROR = 1;

ODEINT_REL_ERROR = 1e-6;

PARTICLE_INITIAL_MASS_FRACTION = 0.01; % dimensionless
% PARTICLE_INITIAL_MASS_FRACTION = 0.025; % dimensionless
% PARTICLE_INITIAL_MASS_FRACTION = 0.05; % dimensionless
% PARTICLE_INITIAL_MASS_FRACTION = 0.075; % dimensionless
%PARTICLE_INITIAL_MASS_FRACTION = 0.1; % dimensionless

PARTICLE_INITIAL_DIAMETER = 85e-9; %m

PARTICLE_SIGMA_G = 1.3;

LOWER_SIZE_SCALE = 0.1; % scaling factor set to initialize smalles particle 
% size in initPBS

DROPLET_DIAMETER = 1.35e-3; %m

DROPLET_TEMPERATURE = 449.15; % K, boiling point of Jet A-1 fuel

DYN_VISC = 0.1978e-3; % dynamic viscosity of n-decane at boiling point, Pa*s

SOLV_DENSITY = 775; %kg/m³, jet fuel A-1

REACTION_RATE = 0; % no precursor present

PARTICLE_DENSITY = 2460; %kg/m³, boron particles with oxide shell

VOLFRAC_CRIT_SURF = 0.64; % 

VOLFRAC_CRIT_CORE = 0.035; % 

NONODES = 80; % number of pivot points for Prakash method

GRID_SPACING_FACTOR = 15; % grid scaling for Prakash method