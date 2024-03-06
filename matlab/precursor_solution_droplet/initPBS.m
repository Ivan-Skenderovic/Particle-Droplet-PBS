%discretize volume domain and initialize particle size distribution
particleMonomerVolume = PARTICLE_MOLAR_MASS/NA/PARTICLE_DENSITY; % m³

% Since the precursor is not included into the coagulation simulation, the
% value of precursorMonomerVolume is not relevant. It is set to particleMonomerVolume
% so that mergeInflowPSD and calcTotalMassPolydisperse can conveniently
% calculate the particle/precursor flux and total particle volume in the simulation
precursorMonomerVolume = particleMonomerVolume; % m³
smallestParticleVolume = particleMonomerVolume; % m³

firstParticleNode = 3;
gridVols = nodalGridVolumes(precursorMonomerVolume, ...
    particleMonomerVolume, ...
    smallestParticleVolume, ...
    NONODES, GRID_SPACING_FACTOR, firstParticleNode);
splitOps = sizeSplittingOperators(gridVols, NONODES); 

smallesParticleDiam = (6*smallestParticleVolume/pi)^(1/3); 
surfaceShellWidth = smallesParticleDiam;

particle_number_concentration = NA*PRECURSOR_CONCENTRATION;

particleConc_surf(1:NONODES) = 0;
particleConc_surf(1) = particle_number_concentration; 

particleConc_core(1:NONODES) = 0;
particleConc_core(1) = particle_number_concentration; 

gridDiams = (6*gridVols/pi).^(1/3);

dropletVolume = pi/6*DROPLET_DIAMETER^3;
totalParticleVolume = particle_number_concentration*smallestParticleVolume*...
    dropletVolume; 
