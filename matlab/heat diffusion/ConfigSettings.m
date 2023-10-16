%%%%%%%%%% Solver Settings %%%%%%%%%%

TIME_START = 1e-10; %s
TIME_END = 1; %s
NOEXPTIMEPOINTS = 100; %

ODEINT_ABS_ERROR = 1e-18; %
ODEINT_REL_ERROR = 1e-18; %
ODEINT_MAX_TIMESCALE_FACTOR = 10; %
ODEINT_MIN_TIMESCALE_FACTOR = 0.0001; %
ODEINT_SAFETY_FACTOR = 0.9; %

NOSTREAMS = 1; %
NOTRAJECTORIES = 1; %
NOSIMS = 1; %
NODROPLETSPERSIM = 1; %
NOPARTICLESPERSIM = 1000; %

%%%%%%%%%% Gas phase %%%%%%%%%%%%

FUEL_OX_MASS_RATIO = 0.208; % stochiometric []
OX_MASS_FRAC = 0.232; % oxidizer mass fraction at droplet surface
PRESSURE_INF = 1; % atm
TEMP_INF = 1700; % K
CARRIER_GAS_VELOCITY = 80; %m/s
COOLING_RATE = 0; %K/s
DILLUTION_FACTOR = 0.05;   % in percent
DILLUTION_TIME = 0.33;	 %s
SAT_START = 1.0001; 
TEMP_START = 1700; %K
TEMP_END = 300; %K

%%%%%%%%%% Finite Difference Configuration %%%%%%%%%%

NOSPATIALGRIDPOINTS = 200; %
MODULE_MAX = 0.3; %
BOUNDARY_CONDITION_TYPE = 3; % 1: insulated 2: temperature dependant heat flux 3: sirignano b.c. for droplet combustion
DROPLET_IDX_SIM_PRINT = 1; %
DROPLET_CUTOFF_DIAMETER = 1e-6; % diameter [m] after instantaneous evaporation is assumed

%%%%%%%%%% Constants %%%%%%%%%%

KB = 1.38064852e-23; %J/K
NA = 6.022140857e23; %1/mol 
R = 8.3144598; %kg*m²/(s²*mol*K)
PI = 3.1415926535;

%%%%%%%%%% Droplets %%%%%%%%%%%%

DROPLET_COAGULATION = 0; % 1 for true
DROPLET_COMBUSTION_MODEL = 3; %1: combustion in quiescent air 2: combustion in convective flow 3: sirignano
DROPLET_RELATIVE_VELOCITY = 80; % m/s	
DROPLET_INITIAL_CONCENTRATION = 1e10; % #/m³
DROPLET_DSD_SHIFT = 0; % [m]
DROPLET_DSD_MASS_SCALING_FACTOR = 0.239261543; %
DROPLET_POP_LOWEST_DIAMETER = 1e-6; %m
DROPLET_POP_MEAN_DIAMETER = 10e-6; %m
DROPLET_POP_HIGHEST_DIAMETER = 1e-3; %m
DROPLET_POP_SIGMA = 1.3; % initially logn distributed
DROPLET_FLAME_TEMPERATURE = 1700; %K
DROPLET_THERMLAY_WIDTHCOEFF = 0.5; %
DROPLET_THERMLAY_HEAT_COND = 0.0398; % W/m/K approximated by air 225°C
DROPLET_INITIAL_TEMPERATURE = 300; %K

%%%%%%%%%%Other Compounds%%%%%

Air_MOLAR_MASS = 0.02897; %kg/mol
Air_DENSITY = 0.348; % kg/m³ at 1000 K 1 bar
Air_HEAT_CAP = 1329.0; % J/kg/K at 1300 K and 1 bar;
Air_HEAT_COND = 4.57e-2; % W/m/K at 600 K and 1 bar 
Air_DYNAMIC_VISCOSITY = 49.0e-6; % Pa*s at 1300 K 1 bar
Air_PRANDTL = 0.7249; %derived from values above for 1300 K and 1 bar