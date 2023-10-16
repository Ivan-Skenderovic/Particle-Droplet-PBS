function [collRateDiff, collRateReg, collRateCirc] = collisionRatesMonodisperse( ...
    mean_particle_diam, kB, temperature, dyn_visc, surface_velocity, shear_rate)

    collRateDiff = 8*kB*temperature/(3*dyn_visc);

    collRateReg = pi.*surface_velocity*mean_particle_diam.^2;
    
    collRateCirc = 4/3*shear_rate*mean_particle_diam.^3;
    
end

