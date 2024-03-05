% Solves precursor thermal decomposition analytically for surface and core
% for simplification purpose, the precursor monomer volume (1st node) is 
% set to the smallest particle volume (3rd node). In case
% smallestParticleVolume is different to particleMonomerVolume volumeRatio
% is not 1 and the concentration needs to be adjusted accordingly.
volumeRatio = gridVols(2)/gridVols(firstParticleNode);

precConcChange_surf = precursorThermalDecomposition(particleConc_surf(1),...
    REACTION_RATE, TIMESTEP);
particleConc_surf(1) = particleConc_surf(1) - precConcChange_surf;
particleConc_surf(firstParticleNode) = particleConc_surf(firstParticleNode) + ...
    volumeRatio*precConcChange_surf;

precConcChange_core = precursorThermalDecomposition(particleConc_core(1),...
    REACTION_RATE, TIMESTEP); 
particleConc_core(1) = particleConc_core(1) - precConcChange_core;
particleConc_core(firstParticleNode) = particleConc_core(firstParticleNode) + ...
    volumeRatio*precConcChange_core;