function [particleConc_surf, particleConc_core, massConservationErrorMerge ] = ...
    mergeInflowPSD(gridVols, r, r_new, shellWidth, shellWidth_new,...
    conc_core, conc_surf)
   
	if (r_new > r)
	%in case of droplet expansion do nothing
	return;
	end
   
    particleVolumeBefore = calcTotalMassPolydisperse(gridVols, conc_core, ...
    conc_surf, shellWidth, r);

    volChange_core = calcShellVolume(r - shellWidth, 0) - ...
        calcShellVolume(r_new - shellWidth_new, 0);
    volFracs_core = conc_core.*gridVols;
    volFluxes = volFracs_core.*volChange_core;
    numberOfParticlesToSurf = volFluxes./gridVols;

    volChange_surf = calcShellVolume(r, r - shellWidth) - ...
        calcShellVolume(r_new, r_new - shellWidth_new);
    volFracs_surf = conc_surf.*gridVols;
    volFluxesToCore = volFracs_surf.*volChange_surf;
    noParticlesToCore = volFluxesToCore./gridVols;

    particleConc_surf(1:length(gridVols)) = conc_surf;
    initialParticleNumber_surf = conc_surf.*...
       calcShellVolume(r, r - shellWidth);
    shellVolume_new = calcShellVolume(r_new, r_new - shellWidth_new);

    particleConc_core(1:length(gridVols)) = conc_core;
    initialParticleNumber_core = conc_core.*...
             calcShellVolume(r - shellWidth, 0);
    coreVolume_new = calcShellVolume(r_new - shellWidth_new, 0);

    for i= 1:length(gridVols) %to noNodes
        
        if shellWidth_new < shellWidth
        %particles flow into core volume        
         particleConc_core(i) = (initialParticleNumber_core(i) + ...
         noParticlesToCore(i))/coreVolume_new;
        else
        %particles flow into surfaceShell volume
         particleConc_surf(i) = (initialParticleNumber_surf(i) + ...
         numberOfParticlesToSurf(i))/shellVolume_new;
        end
        
    end
    
    particleVolumeAfter = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
    particleConc_surf, shellWidth_new, r_new);
 
    massConservationErrorMerge = abs(particleVolumeBefore - particleVolumeAfter)/...
        particleVolumeBefore;
  
end

