disp('Simulation finished. Storing results.');

d_d0 = 2*dropletRadius_new/DROPLET_DIAMETER;       
shellWidth = surfaceShellWidth_new*1e9;
gridDiams_gasPhase = (6.*gridVols_gasPhase(firstParticleNode:end)/pi).^(1/3);
particleConcentrations = particleConc_gas_final(end, firstParticleNode:end);
meanDiam = weightedGeoMean(gridDiams_gasPhase, particleConcentrations);
geostd = sigmag(gridDiams_gasPhase, particleConcentrations);
totalParticleConcentration = sum(particleConcentrations);

% setup results matrix
% export particle number concentration N
output(1, :) = {'d_break/d_0', 'surfaceShellWidth [nm]', 'time [s]', ...
    'mean diam [nm]', 'geo. std [-]', 'N_total [#/m³]'};
output(2, :) = {d_d0, shellWidth, simTime, meanDiam, geostd, ...
    totalParticleConcentration};
output(3, 1:8) = {'particle size', ...
   'N in droplet core at breakup', ...
   'N in droplet surface shell at breakup', ...
   'particle size gas phase', ...
   'N after breakup', ...
   'N final', ...
   'binMidPoints [m]', ...
   'dN/(dlogDp*N_total) [-]'};
output(4, :) = {'[m³]', '[#/m³]', '[#/m³]', '[m³]', '[#/m³]', '[#/m³]',...
    '[m]','[-]'};
output(5 : NONODES + 4, 1) = num2cell(gridVols);
output(5 : NONODES + 4, 2) = num2cell(particleConc_core);
output(5 : NONODES + 4, 3) = num2cell(particleConc_surf);
output(5 : NONODES + 4, 4) = num2cell(gridVols_gasPhase);
output(5 : NONODES + 4, 5) = num2cell(particleConc_gas);
output(5 : NONODES + 4, 6) = num2cell(particleConc_gas_final(end, :));

%setup bin discretization
noBins = 60;
binLowerBound = log10(1e-10);
binUpperBound = log10(1e-4);
binDiams = (6*gridVols_gasPhase(firstParticleNode:end)/pi).^(1/3);

[binMidPoints, binWeights, binWidth] = normalize_psd(...
    binDiams, ...
    particleConcentrations, ...
    noBins, binLowerBound, binUpperBound);

binWeightsNormalized = binWeights./totalParticleConcentration;

output(5 : noBins + 4, 7) = num2cell(binMidPoints);
output(5 : noBins + 4, 8) = num2cell(binWeightsNormalized);

%uncomment to plot a figure
binMidPoints_nm = binMidPoints*1e9;
figure
hold on
box on
plot(binMidPoints_nm, binWeightsNormalized, ...
   'Marker', 'o', 'LineStyle','-.','Color','r','MarkerFaceColor','r');
xlabel('Particle diameter D_P [nm]')
ylabel('dN/(dlogd_P/N_{total}) [-]')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim([1 1e4])
ylim([1e-2 1.5])

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