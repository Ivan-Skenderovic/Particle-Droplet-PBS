
NA = 6.02214086e23; % 1/mol, Avogadro constant

KB = 1.380649e-23; % J/K, Boltzmann constant

TIMESTEP = 1e-7; % s

TIME_END = 1e-3; % s

NOSPATIALGRIDPOINTS = 100; % no of discrete space nodes in finite difference
% calculations

MASS_CONSERVATION_ERROR = 1e-10; % very large particles might "flow" out of 
% the volume grid set by the pivot method

ODEINT_ABS_ERROR_HEATING = 1e-6;

ODEINT_REL_ERROR_HEATING = 1e-3;

ODEINT_ABS_ERROR_COAGULATION = 1e-2;

ODEINT_REL_ERROR_COAGULATION = 1e-6;

PRECURSOR_CONCENTRATION = 500; % mol/m³

DROPLET_DIAMETER = 9e-6; %m

DROPLET_INITIAL_TEMPERATURE = 300; % K

DROPLET_CONCENTRATION = 1e10; % #/m³,  the spray is approximated by a 
% monodisperse droplet population

REACTION_RATE = 1e10; % 1/s

REACTION_RATE_FLAME = 1e10; % 1/s

TEMPERATURE_FLAME = 2300; % K

PARTICLE_MOLAR_MASS = 0.10687; % kg/mol, fe(OH)3

PARTICLE_DENSITY = 3320; % kg/m³, fe(OH)3, bernalite

PARTICLE_DENSITY_GAS_PHASE = 4860; % kg/m³, fe2O3, maghemite

TEMPERATURE_CRIT = 417.15; % K, boiling point of decomposition product

VOLFRAC_CRIT_SURF = 0.1625; % Peclet number is >> 1 for droplet surface

VOLFRAC_CRIT_CORE = 0.035; % Peclet number is 0 for droplet core

NONODES = 320; % number of nodes for Prakash method, set to 80 if 
% computation time is of concenrn

GRID_SPACING_FACTOR = 15; % grid scaling for Prakash method