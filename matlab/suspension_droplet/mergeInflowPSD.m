function [particleConc_surf, massConservationErrorMerge ] = ...
    mergeInflowPSD(gridVols, r, r_new, shellWidth, shellWidth_new, conc_core, conc_surf)
   
    particleVolumeBefore = calcTotalMassPolydisperse(gridVols, conc_core, ...
    conc_surf, shellWidth, r);

    volFracs_core = conc_core.*gridVols;
   
    volChange_core = calcShellVolume(r - shellWidth, 0) - calcShellVolume(r_new - shellWidth_new, 0);

    volFluxes = volFracs_core.*volChange_core;

    numberOfParticlesIncoming = volFluxes./gridVols;

    initialParticleNumber = conc_surf.*calcShellVolume(r, r - shellWidth);

    particleConc_surf(1:length(gridVols)) = 0;

    for i=1:length(gridVols) %to noNodes
         particleConc_surf(i) = (initialParticleNumber(i) + numberOfParticlesIncoming(i))/calcShellVolume(r_new, r_new - shellWidth_new);       
    end
    
     particleVolumeAfter = calcTotalMassPolydisperse(gridVols, conc_core, ...
     particleConc_surf, shellWidth_new, r_new);
 
    massConservationErrorMerge = abs(particleVolumeBefore - particleVolumeAfter)/particleVolumeBefore;
  
end

