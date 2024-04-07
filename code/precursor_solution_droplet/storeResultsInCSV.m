disp('Simulation finished. Storing results.');

d_d0 = 2*dropletRadius_new/DROPLET_DIAMETER;       
shellWidth = surfaceShellWidth_new*1e9;
gridDiams_gasPhase = (6.*gridVols_gasPhase(firstParticleNode:end)/pi).^(1/3);
particleConcentrations = particleConc_gas_final(end, firstParticleNode:end);
meanDiam = lognMedianDiam(gridDiams_gasPhase, particleConcentrations);
geostd = sigmag(gridDiams_gasPhase, particleConcentrations);
totalParticleConcentration = sum(particleConcentrations);

% setup results matrix
output(1, :) = {'d_break/d_0', 'surfaceShellWidth [nm]', 'time [s]', ...
    'mean diam [nm]', 'geo. std [-]', 'total number concentration [#/m³]'};
output(2, :) = {d_d0, shellWidth, simTime, meanDiam, geostd, ...
    totalParticleConcentration};
output(3, 1:6) = {'particle size', ...
   'particle number concentration in droplet core at breakup', ...
   'particle number concentration in droplet surface shell at breakup', ...
   'particle size gas phase', ...
   'particle number concentration after breakup', ...
   'particle number concentration final'};
output(4, :) = {'[m³]', '[#/m³]', '[#/m³]', '[m³]', '[#/m³]', '[#/m³]'};
output(5 : NONODES + 4, 1) = num2cell(gridVols);
output(5 : NONODES + 4, 2) = num2cell(particleConc_core);
output(5 : NONODES + 4, 3) = num2cell(particleConc_surf);
output(5 : NONODES + 4, 4) = num2cell(gridVols_gasPhase);
output(5 : NONODES + 4, 5) = num2cell(particleConc_gas);
output(5 : NONODES + 4, 6) = num2cell(particleConc_gas_final(end, :));

% calc dN/dlogD_p:

outputToTable = cell2table(output);

% set results file name
results_filename = ...
[num2str(PARTICLE_DENSITY_GAS_PHASE,'%2.3g'),'kgm3_',...
num2str(PRECURSOR_CONCENTRATION,'%d'),'M_',...
num2str(DROPLET_DIAMETER*1e6,'%2.3g'),'mu_',...
num2str(TIMESTEP*1e6,'%2.1g'),'musec_',...
num2str(TEMPERATURE_FLAME),'K.csv',]; 

% store as csv file
cd ../../results/
writetable(outputToTable, results_filename); 