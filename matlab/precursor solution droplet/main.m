clear;
clc;
format long e;
close all;
fclose all;

run("ConfigSettings_iron_standard_experiment.m");

% materialproperties
LiqData = LiquidPhaseProperties('2EHA');
GasData = GasPhaseProperties;

initialTemperaturesAndRsquare(1:NOSPATIALGRIDPOINTS) = ...
    DROPLET_INITIAL_TEMPERATURE;
initialTemperaturesAndRsquare(NOSPATIALGRIDPOINTS + 1) = ...
    (DROPLET_DIAMETER*0.5)^2;

options_heatDiffusion = odeset('relTol', ODEINT_REL_ERROR_HEATING, ...
    'absTol', ODEINT_ABS_ERROR_HEATING);
options_coagulation = odeset('relTol', ODEINT_REL_ERROR_COAGULATION, ...
    'absTol', ODEINT_ABS_ERROR_COAGULATION);

initPBS

%% solve PBE for liquid phase until droplet breakup %%%%%%%%%%%%%%%%%%%%%%%
simTime = 0;
isBreaking = false;
noParticleLayers = 1;
isAtLockingPoint = false;

while (~isBreaking) 

simTime = simTime + TIMESTEP;  
dropletRadius = sqrt(initialTemperaturesAndRsquare(end));

% calculate droplet heating and size change 
 odesystemhandle = @(t,Y) HeatDiffusion_1D(t, Y, LiqData, GasData);
 [times, TandRsquare] = ode15s(odesystemhandle, [0 TIMESTEP], ...
     initialTemperaturesAndRsquare, options_heatDiffusion);
 
 dropletRadius_new = sqrt(TandRsquare(end, end));

% particle influx due to droplet radius change 
[particleConc_surf, massConsErrorAdvection] = ...
mergeInflowPSD(gridVols, dropletRadius, dropletRadius_new, ...
surfaceShellWidth, surfaceShellWidth, particleConc_core, particleConc_surf);

% solve PBE for droplet core and surface
coagulationStep

surfaceShellWidth_new = particle_diam_surf_scale*noParticleLayers; 

    if surfaceShellWidth_new > dropletRadius 
    % shellWidth is limited by the droplet radius
        surfaceShellWidth_new = dropletRadius;
    end 

% particle influx due to surfaceShell growth
[particleConc_surf_new, massConsErrorMerge] = ...
mergeInflowPSD(gridVols, dropletRadius_new, dropletRadius_new, surfaceShellWidth, ...
surfaceShellWidth_new, particleConc_core, particleConc_surf);

% check mass conservation
finalParticleVolume = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
particleConc_surf_new, surfaceShellWidth_new, dropletRadius_new); 
massConservationError = abs(finalParticleVolume - totalParticleVolume)/...
    totalParticleVolume;

    if massConservationError > MASS_CONSERVATION_ERROR
        error('Mass conservation failed.')
    end 

%check for droplet events 
isSuperheated = T_core >= TEMPERATURE_CRIT;
volFrac_surf = sum(gridVols(firstParticleNode:end).*...
    particleConc_surf(firstParticleNode:end));

    if (volFrac_surf >= VOLFRAC_CRIT_SURF)
        isAtLockingPoint = true;
        noParticleLayers = noParticleLayers + 1;
    end

% choose one of the following breakup conditions:
 isBreaking = isSuperheated;
% isBreaking = isAtLockingPoint;
% isBreaking = isSuperheated & isAtLockingPoint;

    if (isBreaking) 
        disp(['Breakup event. Particles are now transferred to the ' ...
            'gas phase.']);
       
        transferParticlesToGasphase

        break
    end

initialTemperaturesAndRsquare = TandRsquare(end, :);
particleConc_surf = particleConc_surf_new;
surfaceShellWidth = surfaceShellWidth_new; 

end
%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% solve PBE for gas phase after droplet breakup %%%%%%%%%%%%%%%%%%%%%%%%%% 
%particles shrink due to chemical reaction
gridVols_gasPhase = gridVols*PARTICLE_DENSITY/PARTICLE_DENSITY_GAS_PHASE;
splitOps_gasPhase = sizeSplittingOperators(gridVols, NONODES);
coagConstDiff_fmr = (3/4/pi)^(1/6)*...
    (6*KB*TEMPERATURE_FLAME/PARTICLE_DENSITY_GAS_PHASE)^(0.5);
collRatesDiffFMR = collisionRatesDiff_fmr(coagConstDiff_fmr, ...
    gridVols_gasPhase, NONODES, firstParticleNode);

remainingSimulationTime = TIME_END - simTime;

if (remainingSimulationTime > 0)

    odehandle_gas = @(t,N) solvePBE(t, N, gridVols, collRatesDiffFMR, ...
        splitOps, REACTION_RATE_FLAME, firstParticleNode);	
    [~, particleConc_gas_final] = ode15s(odehandle_gas, ...
        [0 remainingSimulationTime], ...
        particleConc_gas, options_coagulation);      

end
%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

storeResultsInCSV
