clear;
clc;
format long e;
fclose all;
close all;

SETTINGS = 'ConfigSettings_Ohja_Boron.m';
%SETTINGS = 'ConfigSettings_Miglani_Titania.m';
run(SETTINGS);

options = odeset('relTol', ODEINT_REL_ERROR, 'absTol', ODEINT_ABS_ERROR);

r =  0.5*DROPLET_DIAMETER;

initPBS

simTime = 0;
isAtLockingPoint = false;

while (~isAtLockingPoint)

    simTime = simTime + TIMESTEP

    %change of droplet radius using d² Law
    r_new = 0.5*sqrt(4*r^2 - K*TIMESTEP);

    %particle accumulation due to surface regression
    [particleConc_surf, massConsErrorAdvection] = ...
         mergeInflowPSD( gridVols, r, r_new, surfaceShellWidth, surfaceShellWidth, ...
    particleConc_core, particleConc_surf );

    coagulationStep

    surfaceShellWidth_new = particle_diam_surf_scale;
    if surfaceShellWidth_new > r
    %shellWidth is limited by the droplet radius
        surfaceShellWidth_new = r;
    end

    %particle accumulation due to grid adaptation
    [particleConc_surf_new, massConsErrorMerge] = ...
         mergeInflowPSD( gridVols, r_new, r_new, surfaceShellWidth, surfaceShellWidth_new, ...
     particleConc_core, particleConc_surf);

    finalParticleVolume = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
    particleConc_surf_new, surfaceShellWidth_new, r_new);

    massConservationError = ...
        abs(finalParticleVolume - totalParticleVolume)/totalParticleVolume;
    if massConservationError > MASS_CONSERVATION_ERROR
        error('Mass conservation failed.')
    end

    %update for next timestep
    r = r_new;
    volFrac_surf = sum(gridVols.*particleConc_surf);
    particleConc_surf = particleConc_surf_new;
    surfaceShellWidth = surfaceShellWidth_new;

    %terminate simulation at locking point and store results
    if volFrac_surf >= VOLFRAC_CRIT_SURF

       %setup results matrix
       d_d0 = 2*r/DROPLET_DIAMETER;
       shellWidth = surfaceShellWidth_new*1e9;
       output(1, :) = {'d_lock/d_0', 'surfaceShellWidth [nm]', 'time [s]'};
       output(2, :) = {d_d0, shellWidth, simTime};
       output(3, :) = {'particle size', ...
           'particle number concentration in droplet core', ...
           'particle number concentration in droplet surface shell '};
       output(4, :) = {'[m³]', '[#/m³]', '[#/m³]'};
       output(5 : NONODES + 4, 1) = num2cell(gridVols);
       output(5 : NONODES + 4, 2) = num2cell(particleConc_core);
       output(5 : NONODES + 4, 3) = num2cell(particleConc_surf);

       outputToTable = cell2table(output);

       %set results file name
       results_filename = ...
       [num2str(PARTICLE_DENSITY,'%2.3g'),'kgm3_',...
       num2str(PARTICLE_INITIAL_MASS_FRACTION,'%2.3fw'),'_',...
       num2str(DROPLET_DIAMETER*1e6,'%2.3g'),'mu_',...
       num2str(TIMESTEP*1e6,'%2.1g'),'musec_',...
       num2str(K),'K.csv',];

       % store as csv file
       cd ../results/
       writetable(outputToTable, results_filename);

       break

    end

end
