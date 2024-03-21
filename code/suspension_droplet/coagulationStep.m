% calculate particle size change due to coagulation

    coagConstDiff = 2*KB*DROPLET_TEMPERATURE/(3*DYN_VISC);
    collRatesDiff = collisionRatesDiff(coagConstDiff, gridVols, NONODES, ...
        firstParticleNode);

    surface_velocity = (r - r_new)/TIMESTEP; % m/s
    coagConstReg = pi*surface_velocity;
    collRatesReg = collisionRatesReg(coagConstReg, gridVols, NONODES, ...
        firstParticleNode);

% surface volume:
    W_ij_surf = stabilityRatioHyd(particleConc_surf, gridVols, VOLFRAC_CRIT_SURF);
    collisionRatesTotal_surf = (collRatesDiff + collRatesReg)./W_ij_surf;

	odehandle_surf = @(t,N) solvePBE(t, N, gridVols, collisionRatesTotal_surf, ...
        splitOps, REACTION_RATE, firstParticleNode);

    [~, particleConc_surf_new] = ode15s(odehandle_surf, [0 TIMESTEP], ...
        particleConc_surf, options);
    particleConc_surf = particleConc_surf_new(end,:);

    if sum( particleConc_surf(firstParticleNode:end)) > 0
      particle_diam_surf_scale = volWeightedMeanDiam( ...
      gridDiams(firstParticleNode:end), ...
      particleConc_surf(firstParticleNode:end));
    else
      particle_diam_surf_scale = (6*particle_volume/pi)^(1/3);
    end

% core volume:
    W_ij_core = stabilityRatioHyd(particleConc_core, gridVols, VOLFRAC_CRIT_CORE);
    collisionRatesTotal_core = (collRatesDiff)./W_ij_core;

  	odehandle_core = @(t,N) solvePBE(t, N, gridVols, collisionRatesTotal_core, ...
        splitOps, REACTION_RATE, firstParticleNode);

        [t, particleConc_core_new] = ode15s(odehandle_core, [0 TIMESTEP], ...
            particleConc_core, options);
    particleConc_core = particleConc_core_new(end,:);

