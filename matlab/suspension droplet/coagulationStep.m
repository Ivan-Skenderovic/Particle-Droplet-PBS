% calculate particle size change due to coagulation   

    coagConstDiff = 2*KB*DROPLET_TEMPERATURE/(3*DYN_VISC); 
    collRatesDiff = collisionRatesDiff(coagConstDiff, gridVols, NONODES);

    surface_velocity = (r - r_new)/TIMESTEP; %m/s, change of droplet radius over time
    coagConstReg = pi*surface_velocity;  
    collRatesReg = collisionRatesReg(coagConstReg, gridVols, NONODES);
        
% surface volume:  
    W_ij_surf = stabilityRatioHyd(particleConc_surf, gridVols, VOLFRAC_CRIT_SURF);
    collisionRatesTotal_surf = (collRatesDiff + collRatesReg)./W_ij_surf;   
	
    % odehandle_surf = @(t,N) nodal_ode(t, N, gridVols, collisionRatesTotal_surf, splitOps, REACTION_RATE);	
	odehandle_surf = @(t,N) solvePBE(t, N, gridVols, collisionRatesTotal_surf, splitOps, REACTION_RATE);	
	
    [~, particleConc_surf_new] = ode15s(odehandle_surf, [0 TIMESTEP], particleConc_surf, options);    
    particleConc_surf = particleConc_surf_new(end,:);

    if sum( particleConc_surf(3:end)) > 0
      particle_diam_surf_scale = volWeightedMeanDiam(gridDiams(3:end), particleConc_surf(3:end));
    else
      particle_diam_surf_scale = (6*particle_volume/pi)^(1/3);
    end   
    
% core volume:   
    W_ij_core = stabilityRatioHyd(particleConc_core, gridVols, VOLFRAC_CRIT_CORE);    
    collisionRatesTotal_core = (collRatesDiff)./W_ij_core;  

    % odehandle_core = @(t,N) nodal_ode(t, N, gridVols, collisionRatesTotal_core, splitOps, REACTION_RATE);  
	odehandle_core = @(t,N) solvePBE(t, N, gridVols, collisionRatesTotal_core, splitOps, REACTION_RATE);	
	
        [t, particleConc_core_new] = ode15s(odehandle_core, [0 TIMESTEP], particleConc_core, options);    
    particleConc_core = particleConc_core_new(end,:);

    particle_diam_core = meandiamg(gridDiams, particleConc_core);        