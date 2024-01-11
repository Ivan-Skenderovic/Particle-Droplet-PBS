
NA = 6.02214086e23; % 1/mol, Avogadro constant

KB = 1.380649e-23; %J/K, Boltzmann constant

NOSTEPS = 1e6;

TIMESTEP = 1e-6; %s

ODEINT_ABS_ERROR = 1;

ODEINT_REL_ERROR = 1e-6;

PRECURSOR_CONCENTRATION = 50; % mol/m³

PARTICLE_INITIAL_MASS_FRACTION = 0.0; % dimensionless

PARTICLE_INITIAL_DIAMETER = 0.0; % m

INIT_LOGN = 0;

DROPLET_DIAMETER = 18e-6; %m

DROPLET_VOLUME_CUTOFF = 1; 

DROPLET_TEMPERATURE = 501.2; % K, boiling point of 2-EHA

K = 2.0e-6;% droplet surface regression rate, m²/s

REACTION_RATE = 1e10;

DYN_VISC = 7.8e-3; % dynamic viscosity of solvent, Pa*s

SOLV_DENSITY = 910; %kg/m³, 2-EHA

RELATIVE_PERMITTIVITY = 58.5; % relative permittivity of solvent (here formic acid as substitute)

SHELL_SCALEFACTOR = 1; % monolayers required to build a stable particle shell

%PARTICLE_MOLAR_MASS = 0.159687; %kg/mol, fe2O3

%PARTICLE_MOLAR_MASS = 0.071844; %kg/mol, feO

%PARTICLE_MOLAR_MASS = 0.08885; %kg/mol, feO-OH

PARTICLE_MOLAR_MASS = 0.08986; %kg/mol, fe(OH)2

%PARTICLE_MOLAR_MASS = 0.10687; %kg/mol, fe(OH)3

%PARTICLE_DENSITY = 4860; %kg/m³, maghemite

%PARTICLE_DENSITY = 5740; %kg/m³, FeO

%PARTICLE_DENSITY = 7874; %kg/m³, Fe

%PARTICLE_DENSITY = 4250; %kg/m³, FeO-OH

%PARTICLE_DENSITY = 3400; %kg/m³, Fe(OH)2

PARTICLE_DENSITY = 3320; %kg/m³,  fe(OH)3

PARTICLE_OXYGEN_STOCHIOMETRY = 1;

VOLFRAC_CRIT = 0; % depends on fractal dimension of particles

BIMODAL_SPHERE_PACKING = 1; %volume fraction achieved by bimodal psd, set to 1 if unused

NONODES = 80; % number of nodes for Prakash method

GRID_SPACING_FACTOR = 15; % grid scaling for Prakash method