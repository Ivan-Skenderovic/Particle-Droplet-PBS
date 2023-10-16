function totalMass = calcTotalMassPolydisperse(gridVols, conc_core, conc_surf, surfaceShellWidth, r) 

    coreVolume = calcShellVolume( r - surfaceShellWidth, 0);
    coreMass = sum(gridVols.*conc_core)*coreVolume;

    surfaceShellVolume = calcShellVolume( r, r - surfaceShellWidth );
    surfaceShellMass = sum(gridVols.*conc_surf)*surfaceShellVolume;

    totalMass = coreMass + surfaceShellMass;

end

