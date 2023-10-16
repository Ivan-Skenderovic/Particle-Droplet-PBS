clear;
clc;
format long e;
fclose all;
close all;

SETTINGSNAME = 'ConfigSettings_Ohja_Boron.m';
run(SETTINGSNAME);

r =  0.5*DROPLET_DIAMETER;
droplet_volume = pi/6*DROPLET_DIAMETER^3; % m³
droplet_mass = droplet_volume*SOLV_DENSITY;

initPBS

s_massConservationCheck(1:NOSTEPS) = 0;
s_surfaceConcentration(1:NOSTEPS) = 0;
s_volFrac(1:NOSTEPS) = 0;
s_volFracCheck(1:NOSTEPS) = 0;
s_surfaceShellWidth(1:NOSTEPS) = 0;
s_radius(1:NOSTEPS) = 0;
s_particle_diam(1:NOSTEPS) = 0;
s_collisionRate(1:NOSTEPS) = 0;
s_collisionRate_diff(1:NOSTEPS) = 0;
s_collisionRate_reg(1:NOSTEPS) = 0;
s_collisionRate_circ(1:NOSTEPS) = 0;
s_DLVO(1:NOSTEPS) = 0;
s_r_avg(1:NOSTEPS) = 0;
s_meandiamg(1:NOSTEPS) = 0;
s_Peclet(1:NOSTEPS) = 0;
s_solidVolumeSaturation(1:NOSTEPS) = 0;
s_massConservationCheck(1:NOSTEPS) = 0;
s_massConservationMergeCheck(1:NOSTEPS) = 0;
s_massConservationAdvectionCheck(1:NOSTEPS) = 0;
s_massConservationSolverCheck(1:NOSTEPS) = 0;
s_concLastNode_surf(1:NONODES) = 0;
s_concLastNode_core(1:NONODES) = 0;

isLocked = false;
lockingPoint = 0;

%breakage model
dropletInitialVolume = pi/6*DROPLET_DIAMETER^3;
particleConcAccumulated(1:NONODES) = 0;
particleConcAdded(1:NONODES) = 0;
particleVolAccumulated = 0;

simTime = 0;

tic

for i = 1:NOSTEPS    

    %time update                      
    simTime = simTime + TIMESTEP;
    
    %droplet heating, calculate surface regression rate K
    K = 0.8e-6;

    %droplet radius update
    r_new = 0.5*sqrt(4*r^2 - K*TIMESTEP);  
   
    %particle accumulation due to surface regression
     [particleConc_surf, massConsErrorAdvection] = ...
         mergeInflowPSD( gridVols, r, r_new, surfaceShellWidth, surfaceShellWidth, ...
     particleConc_core, particleConc_surf);
       
    coagulationStep
    
    %update surface cell size
    surfaceShellWidth_new = particle_diam_surf_scale; 
    if surfaceShellWidth_new > r
        surfaceShellWidth_new = r;
    end 
    
    %particle accumulation due to grid adaptation
    [particleConc_surf_new, massConsErrorMerge] = ...
         mergeInflowPSD( gridVols, r_new, r_new, surfaceShellWidth, surfaceShellWidth_new, ...
     particleConc_core, particleConc_surf);
   
     %mass conservation check
     finalParticleVolume = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
     particleConc_surf_new, surfaceShellWidth_new, r_new);
     
    %update for next timestep
    r = r_new;    
    volFrac_surf = sum(gridVols.*particleConc_surf);
    particleConc_surf = particleConc_surf_new;
    surfaceShellWidth = surfaceShellWidth_new;  
    
    %particle shell formation check
     if volFrac_surf >= VOLFRAC_CRIT
        if isLocked == false
            
           workspace_vars = ...
          ['Ohja_',num2str(PARTICLE_INITIAL_MASS_FRACTION,'%2.1fw'),'_',...
           num2str(DROPLET_DIAMETER*1e6,'%2.1g'),'mu_',...
           num2str(TIMESTEP*1e6,'%2.1g'),'musec_',...
           num2str(K),'K.mat',];
           save(workspace_vars);  

           lockingDiameter = 2*r;               
           dropletLifeTimeCheck = pi/6*lockingDiameter^3/dropletInitialVolume;
            
            if dropletLifeTimeCheck > DROPLET_VOLUME_CUTOFF
                %break droplet and accumulate particles
                shellVolumeAtBreakup = calcShellVolume(r, r - surfaceShellWidth);
                
                particleConcAdded = shellVolumeAtBreakup*particleConc_surf;
 
                particleConcAccumulated = particleConcAccumulated + particleConcAdded;
                
                %set size of child droplet
                r = r - surfaceShellWidth;
                %set particle concentrations in child droplet shell
                particleConc_surf = particleConc_core;
                
                %recalculate total droplet mass for statistics
                 finalParticleVolume = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
                     particleConc_surf, surfaceShellWidth, r);
                
            else
                simTime
                
                toc
                
                coreVolume = calcShellVolume(r - surfaceShellWidth, 0);
                surfaceShellVolume = calcShellVolume(r, r - surfaceShellWidth);
                
                particleConcAdded = coreVolume*particleConc_core + surfaceShellVolume*particleConc_surf;                
                particleConcAccumulated = particleConcAccumulated + particleConcAdded;                
                
                lockingPoint = i;
                isLocked = true;
                K = 0;                
                    
                %plotting
               
                finalDiameter = 2*r_new                               
                
                break
            end
            
            particleVolAccumulated = particleVolAccumulated + sum(gridVols.*particleConcAdded);
            particleConcAdded = 0;
            
        end
    end
              
    massConservationError = abs((finalParticleVolume + particleVolAccumulated) - totalParticleVolume)/totalParticleVolume;
    if massConservationError > 1e-8
        error('Mass conservation failed.')
    end    
    
    statistics_pbe
   
end
