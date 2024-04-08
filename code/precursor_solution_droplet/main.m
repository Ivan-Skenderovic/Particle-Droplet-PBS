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

        % calculate total iron mole amount after breakup and before
        % decomposition reaction
        volPerMole_feOH3 = PARTICLE_MOLAR_MASS/PARTICLE_DENSITY; % m³/mole
        moleAmount_fe = totalParticleVolumeAfterBreakup/volPerMole_feOH3;

        break
    end

initialTemperaturesAndRsquare = TandRsquare(end, :);
particleConc_surf = particleConc_surf_new;
surfaceShellWidth = surfaceShellWidth_new; 

end
%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% solve PBE for gas phase after droplet breakup and account for particle % 
%% volume change due to chemical reactions: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) precursor decomposition reaction in gas phase:
% 2*Fe-(eha)3 + [...] -> Fe2O3 + [...] 
% mole amount is halved:
stoichiometry_factor = 2;
particleConc_gas(1) = particleConc_gas(1)/stoichiometry_factor; 
% 2) particle thermal decomposition in gas phase:
% 2*Fe(OH)3 -> Fe2O3 + 3*H2O; H2O evaporates from particle

M_maghemite = 159.69e-3; % kg / mole, Molar Mass
volPerMole_maghemite = M_maghemite/PARTICLE_DENSITY_GAS_PHASE; % m³/mole

reactionVolumeRatio = stoichiometry_factor*volPerMole_feOH3/...
    volPerMole_maghemite;

gridVols_gasPhase = gridVols/reactionVolumeRatio;
monomerVolume_gasPhase = volPerMole_maghemite/NA;
gridVols_gasPhase(1:2) = monomerVolume_gasPhase;
% Prevent inclusion of particles smaller than monomer volume. In
% the solvePBE function this is done using the nucleation split
for i = firstParticleNode:length(gridVols_gasPhase)
    if gridVols_gasPhase(i) < monomerVolume_gasPhase
        newVolumeRatio = gridVols_gasPhase(i)/monomerVolume_gasPhase;
        particleConc_gas(i) = particleConc_gas(i)*newVolumeRatio;
        gridVols_gasPhase(i) = monomerVolume_gasPhase;
    end
end

splitOps_gasPhase = sizeSplittingOperators(gridVols_gasPhase, NONODES);
coagConstDiff_fmr = (3/4/pi)^(1/6)*...
    (6*KB*TEMPERATURE_FLAME/PARTICLE_DENSITY_GAS_PHASE)^(0.5);
collRatesDiffFMR = collisionRatesDiff_fmr(coagConstDiff_fmr, ...
    gridVols_gasPhase, NONODES, firstParticleNode);

remainingSimulationTime = TIME_END - simTime;
if (remainingSimulationTime > 0)

    odehandle_gas = @(t,N) solvePBE(t, N, gridVols_gasPhase, collRatesDiffFMR, ...
        splitOps_gasPhase, REACTION_RATE_FLAME, firstParticleNode);	

    [~, particleConc_gas_final] = ode15s(odehandle_gas, ...
        [0 remainingSimulationTime], ...
        particleConc_gas, options_coagulation);      

end

% Mole amount of Fe is conserved in the reaction, particle mass and volume
% are not. Mole amount conservation and mass loss are checked to confirm
% the correctness of the calculation.
particleVolume_gasPhase = sum(particleConc_gas_final(end,:).*gridVols_gasPhase);
moleAmount_fe_gasPhase = particleVolume_gasPhase/...
    volPerMole_maghemite*stoichiometry_factor;
% checkMoleAmount needs to be 1.
checkMoleAmount = moleAmount_fe/moleAmount_fe_gasPhase;

% check mass loss from evaporation of H2O
particleMass_feOH3 = totalParticleVolumeAfterBreakup*PARTICLE_DENSITY;
stoichiometry_factor_h2O = 3/2;
M_h2o = 18.01528e-3; % kg / mole
totalMass_h2O = moleAmount_fe*stoichiometry_factor_h2O*M_h2o;
particleMass_maghemite = particleVolume_gasPhase*...
    PARTICLE_DENSITY_GAS_PHASE;
% checkMassLoss needs to be 1. Rounding errors from Molar masses
% are expected.
checkMassLoss = (particleMass_feOH3 - totalMass_h2O)/...
    particleMass_maghemite; 

%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

storeResultsInCSV
