% calculate particle size distribution evolution due to coagulation 

% surface volume:
    T_surf = TandRsquare(end, end-1);
    
    coagConstDiffSurface = 2*KB*T_surf/(3*LiqData.dyn_visc); 
    collRatesDiffSurface = collisionRatesDiff(coagConstDiffSurface, gridVols, NONODES);
    surface_velocity = (dropletRadius - dropletRadius_new)/TIMESTEP; 
    coagConstReg = pi*surface_velocity;  
    collRatesReg = collisionRatesReg(coagConstReg, gridVols, NONODES);

    W_ij_surf = stabilityRatioHyd(particleConc_surf, gridVols, VOLFRAC_CRIT_SURF);
    collisionRatesTotal_surf = (collRatesDiffSurface + collRatesReg)./W_ij_surf;   
	
    odehandle_surf = @(t,N) solvePBE(t, N, gridVols, collisionRatesTotal_surf, ...
        splitOps, REACTION_RATE, firstParticleNode);	
    [~, particleConc_surf_new] = ode15s(odehandle_surf, [0 TIMESTEP], ...
        particleConc_surf, options_coagulation);    
    particleConc_surf = particleConc_surf_new(end,:);

    if sum( particleConc_surf(firstParticleNode:end)) > 0
      particle_diam_surf_scale = volWeightedMeanDiam(...
          gridDiams(firstParticleNode:end), particleConc_surf(firstParticleNode:end) );
    else
      particle_diam_surf_scale = (6*smallestParticleVolume/pi)^(1/3);
    end   
    
% core volume:   
    Temperatures = TandRsquare(end, 1:NOSPATIALGRIDPOINTS);
    DropletDiameter_new = 2*dropletRadius_new;
    T_core = volWeightedMeanTemperature(Temperatures, DropletDiameter_new);

    coagConstDiffCore = 2*KB*T_core/(3*LiqData.dyn_visc); 
    collRatesDiffCore = collisionRatesDiff(coagConstDiffCore, gridVols, NONODES);

    W_ij_core = stabilityRatioHyd(particleConc_core, gridVols, VOLFRAC_CRIT_CORE);    
    collisionRatesTotal_core = (collRatesDiffCore)./W_ij_core;  

    odehandle_core = @(t,N) solvePBE(t, N, gridVols, collisionRatesTotal_core,...
    splitOps, REACTION_RATE, firstParticleNode);   
    [~, particleConc_core_new] = ode15s(odehandle_core, [0 TIMESTEP],...
        particleConc_core, options_coagulation);    
    particleConc_core = particleConc_core_new(end,:);       