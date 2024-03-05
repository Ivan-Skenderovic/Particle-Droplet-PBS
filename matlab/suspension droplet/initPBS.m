%discretize volume domain and initialize particle size distribution

% if PARTICLE_INITIAL_MASS_FRACTION > 0
    particle_mass_total = droplet_mass*PARTICLE_INITIAL_MASS_FRACTION;
    particle_diam = PARTICLE_INITIAL_DIAMETER;
    particle_volume = pi/6*particle_diam^3;
    particle_mass = particle_volume*PARTICLE_DENSITY;
    particle_number = particle_mass_total/particle_mass;    
    particle_number_concentration = particle_number/droplet_volume;
% else
%     %particles are formed from precursor molecules
%     particle_volume = PARTICLE_MOLAR_MASS/NA/PARTICLE_DENSITY; %m, initially a monomer
%     particle_diam = (6*particle_volume/pi)^(1/3);
%     particle_number_concentration = NA*PRECURSOR_CONCENTRATION/PARTICLE_OXYGEN_STOCHIOMETRY;
% end   

particle_volume_fraction = particle_number_concentration*particle_volume;
totalParticleVolume = particle_volume_fraction*droplet_volume;
surfaceShellWidth = particle_diam;

volFrac_core = particle_volume_fraction;
volFrac_surf = particle_volume_fraction;

%% setup nodal method

%if INIT_LOGN %initial PSD ist log normal 
    precursorMonomerVolume = pi/6*(PARTICLE_INITIAL_DIAMETER*0.1)^3;
    particleMonomerVolume = pi/6*(PARTICLE_INITIAL_DIAMETER*0.1)^3;
    smallestParticleVolume = pi/6*(PARTICLE_INITIAL_DIAMETER*0.1)^3; 
    
    gridVols = nodalGridVolumes(precursorMonomerVolume, particleMonomerVolume, smallestParticleVolume, ...
        NONODES, GRID_SPACING_FACTOR);
    gridDiams = (6*gridVols/pi).^(1/3);
    
   sigma_g = 1.3; 
   particle_volume_total = particle_mass_total/PARTICLE_DENSITY;     
   particleConcentrations = createLognPsd(NONODES, gridDiams, PARTICLE_INITIAL_DIAMETER, sigma_g, particle_volume_total, ...
   droplet_volume, 1e-5);
	
   %scale to set correct particle volume
   particle_volume_initLogn = sum(particleConcentrations.*gridVols)*droplet_volume;
   scaling = particle_volume_total/particle_volume_initLogn;
   particleConcentrations = particleConcentrations*scaling;

   volumeCheck = abs(sum(particleConcentrations.*gridVols)*droplet_volume - particle_volume_total)/particle_volume_total;
   if volumeCheck > 1e-15
       error('particle volume init wrong!')
   end
   
   particleConc_core = particleConcentrations;
   particleConc_surf = particleConcentrations;
   
% else %initial PSD is monodisperse
%     precursorMonomerVolume = particle_volume;
%     particleMonomerVolume = particle_volume;
%     smallestParticleVolume = particle_volume;
%     particle_diam_surf = PARTICLE_INITIAL_DIAMETER;
% 
%     gridVols = nodalGridVolumes(precursorMonomerVolume, particleMonomerVolume, smallestParticleVolume, ...
%         NONODES, GRID_SPACING_FACTOR);
%     gridDiams = (6*gridVols/pi).^(1/3);
%     
%     particleConc_surf(1:NONODES) = 0;
%     particleConc_surf(1) = particle_number_concentration; 
% 
%     particleConc_core(1:NONODES) = 0;
%     particleConc_core(1) = particle_number_concentration; 
% end

splitOps = sizeSplittingOperators(gridVols, NONODES); 
particle_diam_core = meandiamg(gridDiams, particleConc_core); 
%%