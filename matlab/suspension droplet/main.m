clear;
clc;
format long e;
fclose all;
close all;

SETTINGS = 'ConfigSettings_Ohja_Boron.m';
%SETTINGS = 'ConfigSettings_Miglani_Titania.m';
run(SETTINGS);

options = odeset('relTol', ODEINT_REL_ERROR, 'absTol', ODEINT_ABS_ERROR);

LiqData = LiquidPhaseProperties('Octane');
GasData = GasPhaseProperties;

r =  0.5*DROPLET_DIAMETER;
droplet_volume = pi/6*DROPLET_DIAMETER^3; % m³
droplet_mass = droplet_volume*SOLV_DENSITY;

initPBS

initialTemperaturesAndR2(1:NO_FDGRIDPOINTS) = 0;
initialTemperaturesAndR2(1:end-1) = 300; %set initial temperature as room temperature
initialTemperaturesAndR2(end) = (0.5*DROPLET_DIAMETER)^2;

simTime = 0;

for i = 1:NOSTEPS    

    %time update                      
    simTime = simTime + TIMESTEP;   

    calcTemperatureAndSizeChange
  
    %particle accumulation due to surface regression
%      [particleConc_surf, massConsErrorAdvection] = ...
%          mergeInflowPSD( gridVols, r, r_new, surfaceShellWidth, surfaceShellWidth, ...
%      particleConc_core, particleConc_surf );
       
    %coagulationStep
    
    %update surface cell size
%     surfaceShellWidth_new = particle_diam_surf_scale; 
%     if surfaceShellWidth_new > r
%         surfaceShellWidth_new = r;
%     end 
%     
%     %particle accumulation due to grid adaptation
%     [particleConc_surf_new, massConsErrorMerge] = ...
%          mergeInflowPSD( gridVols, r_new, r_new, surfaceShellWidth, surfaceShellWidth_new, ...
%      particleConc_core, particleConc_surf);
%    
%     %calculate for mass conservation check
%     finalParticleVolume = calcTotalMassPolydisperse(gridVols, particleConc_core, ...
%     particleConc_surf_new, surfaceShellWidth_new, r_new);
 
    %update for next timestep
    r = r_new;    
%     volFrac_surf = sum(gridVols.*particleConc_surf);
%     particleConc_surf = particleConc_surf_new;
%     surfaceShellWidth = surfaceShellWidth_new;  
%     
    %terminate simulation at locking point and store results
%     if volFrac_surf >= VOLFRAC_CRIT 
%             
%        workspace_vars = ...
%        ['Ohja_',num2str(PARTICLE_INITIAL_MASS_FRACTION,'%2.3fw'),'_',...
%        num2str(DROPLET_DIAMETER*1e6,'%2.3g'),'mu_',...
%        num2str(TIMESTEP*1e6,'%2.1g'),'musec_',...
%        num2str(K),'K.mat',];
%        save(workspace_vars);  
%        
%        break
% 
%     end
%               
%     massConservationError = ...
%         abs(finalParticleVolume - totalParticleVolume)/totalParticleVolume;
%     if massConservationError > 1e-10
%         error('Mass conservation failed.')
%     end    
   
end
