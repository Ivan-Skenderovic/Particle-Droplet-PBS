% calculates new particle concentrations in the gas phase after breakup
% event

dropletRadius_inner = dropletRadius_new - surfaceShellWidth_new;
dropletShellVol = calcShellVolume(dropletRadius_new, ...
    dropletRadius_inner);
DropletDiameter_new = 2*dropletRadius_new;
dropletVol = pi/6*DropletDiameter_new^3;
dropletCoreVol = dropletVol - dropletShellVol;

    if (isAtLockingPoint)
    % Particles in the surface shell enter the gas phase as a single ...
    % agglomerate. Precursor enters individually.
    
        % add precursor and particles from core
        particleConc_gas = (dropletCoreVol.*particleConc_core)*...
            DROPLET_CONCENTRATION;
    
        % add surface agglomerate by splitting between adjacent nodes
        surfAgglomerateVol = sum( gridVols(firstParticleNode:end).*...
            particleConc_surf_new(firstParticleNode:end)).*dropletShellVol;
        surfAgglomerateConc = DROPLET_CONCENTRATION;
    
        rightIdx = find(gridVols > surfAgglomerateVol, 1, 'first');
    
        volumeSplit = (surfAgglomerateVol - gridVols(rightIdx-1))/...
            (gridVols(rightIdx) - gridVols(rightIdx-1));
    
        particleConc_gas(rightIdx) =  particleConc_gas(rightIdx) + ...
            volumeSplit*surfAgglomerateConc;
        particleConc_gas(rightIdx - 1) = particleConc_gas(rightIdx - 1) + ...
            (1-volumeSplit)*surfAgglomerateConc;

        % add precursor from surface shell
        surfPrecursorConc = particleConc_surf_new(1)*dropletShellVol*...
            DROPLET_CONCENTRATION;    
        particleConc_gas(1) = particleConc_gas(1) + surfPrecursorConc;
    
    else
    % Particles and precursor enter the gas phase individually.
        particleConc_gas = (dropletCoreVol.*particleConc_core + ...
        dropletShellVol.*particleConc_surf_new)*DROPLET_CONCENTRATION;
    end

%check mass conservation
finalParticleVolume_new = finalParticleVolume*DROPLET_CONCENTRATION;
totalParticleVolumeAfterBreakup = sum(gridVols.*particleConc_gas);

massConservationError = abs(finalParticleVolume_new - ...
    totalParticleVolumeAfterBreakup)/totalParticleVolumeAfterBreakup;

    if massConservationError > MASS_CONSERVATION_ERROR
        error('Mass conservation failed.')
    end 