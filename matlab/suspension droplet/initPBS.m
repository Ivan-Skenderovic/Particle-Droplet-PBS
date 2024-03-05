% discretize volume domain and initialize particle size distribution

particle_mass_total = droplet_mass*PARTICLE_INITIAL_MASS_FRACTION;
particle_diam = PARTICLE_INITIAL_DIAMETER;
particle_volume = pi/6*particle_diam^3;
particle_mass = particle_volume*PARTICLE_DENSITY;
particle_number = particle_mass_total/particle_mass;    
particle_number_concentration = particle_number/droplet_volume; 

particle_volume_fraction = particle_number_concentration*particle_volume;
totalParticleVolume = particle_volume_fraction*droplet_volume;
surfaceShellWidth = particle_diam;

volFrac_core = particle_volume_fraction;
volFrac_surf = particle_volume_fraction;

% set smallest particle volume, in suspension droplet no precursor or particle monomers are present
% they are included here in case the user wants to simulate multi component precursor/particle systems 
precursorMonomerVolume = pi/6*(PARTICLE_INITIAL_DIAMETER*LOWER_SIZE_SCALE)^3;
particleMonomerVolume = pi/6*(PARTICLE_INITIAL_DIAMETER*LOWER_SIZE_SCALE)^3;
smallestParticleVolume = pi/6*(PARTICLE_INITIAL_DIAMETER*LOWER_SIZE_SCALE)^3; 

gridVols = nodalGridVolumes(precursorMonomerVolume, particleMonomerVolume, smallestParticleVolume, ...
	NONODES, GRID_SPACING_FACTOR);
gridDiams = (6*gridVols/pi).^(1/3);

particle_volume_total = particle_mass_total/PARTICLE_DENSITY;     
particleConcentrations = createLognPsd(NONODES, gridDiams, PARTICLE_INITIAL_DIAMETER, PARTICLE_SIGMA_G, particle_volume_total, ...
droplet_volume, 1e-5);

%s cale to set correct particle volume
particle_volume_initLogn = sum(particleConcentrations.*gridVols)*droplet_volume;
scaling = particle_volume_total/particle_volume_initLogn;
particleConcentrations = particleConcentrations*scaling;

	volumeCheck = abs(sum(particleConcentrations.*gridVols)*droplet_volume - particle_volume_total)/particle_volume_total;
	if volumeCheck > 1e-15
	   error('particle volume initialization went wrong!')
	end

particleConc_core = particleConcentrations;
particleConc_surf = particleConcentrations;
   

splitOps = sizeSplittingOperators(gridVols, NONODES); 
particle_diam_core = meandiamg(gridDiams, particleConc_core); 
