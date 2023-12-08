% calculate particle size change due to coagulation   

particleVolumeSolverBefore = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
     particleConc_surf, surfaceShellWidth, r_new);
      
    coagConstDiff = 2*KB*DROPLET_TEMPERATURE/(3*DYN_VISC); % 1/m³/s, collision rate due to Brownian diffusion
    collRatesDiff = collisionRatesDiff(coagConstDiff, gridVols, NONODES);

    surface_velocity = (r - r_new)/TIMESTEP; %m/s, dRdt
    coagConstReg = pi*surface_velocity; % 
    collRatesReg = collisionRatesReg(coagConstReg, gridVols, NONODES);
        
% surface volume:  
    W_ij_surf = stabilityRatioHyd(particleConc_surf, gridDiams, VOLFRAC_CRIT);
    collisionRatesTotal_surf = (collRatesDiff + collRatesReg)./W_ij_surf;   
    odehandle_surf = @(t,N) nodal_ode(t, N, gridVols, collisionRatesTotal_surf, splitOps, REACTION_RATE);
    options_surf = odeset('absTol', ODEINT_ABS_ERROR, 'relTol', ...
    ODEINT_REL_ERROR, 'Nonnegative', 1);
    
    [t, particleConc_surf_new] = ode15s(odehandle_surf, [0 TIMESTEP], particleConc_surf, options_surf);    
    particleConc_surf = particleConc_surf_new(end,:);

    if sum( particleConc_surf(2:end)) > 0
      particle_diam_surf_scale = volWeightedMeanDiam(gridDiams(3:end), particleConc_surf(3:end));
    else
      particle_diam_surf_scale = (6*particle_volume/pi)^(1/3);
    end   
    
% core volume:   
    W_ij_core = stabilityRatioHyd(particleConc_core, gridDiams, 1/19);    
    collisionRatesTotal_core = (collRatesDiff)./W_ij_core;  

    odehandle_core = @(t,N) nodal_ode(t, N, gridVols, collisionRatesTotal_core, splitOps, REACTION_RATE);
    options_core = odeset('absTol', ODEINT_ABS_ERROR, 'relTol', ...
    ODEINT_REL_ERROR, 'Nonnegative', 1);
   
        [t, particleConc_core_new] = ode15s(odehandle_core, [0 TIMESTEP], particleConc_core, options_core);    
    particleConc_core = particleConc_core_new(end,:);

    particle_diam_core = meandiamg(gridDiams, particleConc_core);         


particleVolumeSolverAfter = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
     particleConc_surf, surfaceShellWidth, r_new);
 
 massConsErrorSolver = abs(particleVolumeSolverBefore - particleVolumeSolverAfter)/particleVolumeSolverBefore;