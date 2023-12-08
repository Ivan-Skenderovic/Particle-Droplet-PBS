%% plots droplet life time related quantities
%load('iron_11microm_050M.mat')

time = (1:NOSTEPS)*TIMESTEP;

figure(1)
subplot(2,2,1)
box on
hold on
plot(time(1:i), s_massConservationCheck(1:i));
plot(time(1:i), s_massConservationMergeCheck(1:i));
plot(time(1:i), s_massConservationAdvectionCheck(1:i));
plot(time(1:i), s_massConservationSolverCheck(1:i));
set(gca,'yscale','log');
ylabel('massConservation error')
xlabel('time [s]')
legend('total','merge step particle growth','merge step advection','solver')

subplot(2,2,2)
box on
hold on
plot(time(1:i), s_solidVolumeSaturation(1:i));
set(gca,'yscale','log');
ylabel('solidVolumeSaturation')
xlabel('time [s]')

subplot(2,2,3)
box on
hold on
plot(time(1:i), s_surfaceShellWidth(1:i));
plot(time(1:i), s_radius(1:i));
set(gca,'yscale','log');
ylabel('surfaceShellWidth [m]')
xlabel('time [s]')
legend('shellWidth','radius')

% subplot(2,2,4)
% box on
% hold on
% plot(time(1:i), s_concLastNode_surf(1:i));
% plot(time(1:i), s_concLastNode_core(1:i));
% set(gca,'yscale','log');
% ylabel('particle concentration last node [1/m³]')
% xlabel('time [s]')
% legend('surface','core')

% experiment_fe_EtOH = load('Li_Fe_EtOH.csv');
% time_experiment_fe_EtOH = experiment_fe_EtOH(:,1);
% size_experiment_fe_EtOH = experiment_fe_EtOH(:,2);

droplet_LifeTime = 1  - ( s_radius./s_radius(1) ).^3;
collisionRateRatio = s_collisionRate_reg./s_collisionRate_diff;
normedVolumeFraction = s_volFrac./s_volFracCheck;

sim_normed_droplet_diameter = (s_radius./s_radius(1)).^2;
active_droplet_diameter =  sim_normed_droplet_diameter( 1:find( sim_normed_droplet_diameter < 1e-6, 1 ) );
size_active = length(active_droplet_diameter);
sim_normed_droplet_diameter(size_active:end) = sim_normed_droplet_diameter(size_active-1);
sim_normed_time = (1:NOSTEPS)*TIMESTEP./((2*s_radius(1)).^2*1e6); % s/m²
droplet_LifeTime = droplet_LifeTime(1:size_active-1);

% figure(4)
% box on;
% hold on;
% % yyaxis left
% % plot(droplet_LifeTime, normedVolumeFraction(1:size_active-1));
% xlabel('1 - d^3/d_0^3')
% % ylabel('\phi / \phi_{crit}')
% % xlim([0 1])
% % yyaxis right
% plot(droplet_LifeTime, s_r_avg(1:size_active-1));
% ylabel(' norm. separation distance 2r/d_g [-]')

% figure(100)
% box on;
% hold on;
% plot(sim_normed_time, s_radius.^2./s_radius(1).^2);
% xlabel('t/d_0^2 [\mum/\mus^2]')
% ylabel('d^2 / d^2_0 [-]')

% figure(1)
% box on;
% hold on;
% yyaxis left
% plot(droplet_LifeTime, s_particle_diam_surf(1:size_active-1)*1e9);
% plot(droplet_LifeTime, s_particle_diam_surf_scale(1:size_active-1)*1e9);
% set(gca,'yscale','log')
% ylabel('average diameter / nm')
% yyaxis right
% plot(droplet_LifeTime, normedVolumeFraction(1:size_active-1));
% ylabel('\phi / \phi_{crit}')
% xlabel('1 - d^3/d_0^3')
% legend('\gamma d_{4,3}','d_{5,0}')

%plot experiment
% experiment_data_stodt = load('EtOH_2EHA_INN_DropletDiameterEvoloution_Stodt.csv');
% %experiment_data = load('Sim_Al_nano_0.06.csv');
% time_experiment_stodt = experiment_data_stodt(:,1); %in ms
% size_experiment_stodt = mean(experiment_data_stodt(:,2:end),2); % in micrometer
% %convert to s/mm²
% time_experiment_stodt = time_experiment_stodt*1e-3/( (size_experiment_stodt(1)*1e-6)^2*1e6);
% %norm by initial size
% size_experiment_stodt = size_experiment_stodt./size_experiment_stodt(1);
% 
% experiment_data_li = load('Li_Fe_2EHA_EtOH.csv');
% time_experiment_li = experiment_data_li(:,1); %in micro s/micro m²
% size_experiment_li = experiment_data_li(:,2); % dimensionless d²
% timeOffset = 0.0;

% figure(10)
% box on
% hold on
% plot( time_experiment_li, size_experiment_li, 'LineStyle','none','Marker', 'o','MarkerSize',2,'Color','k','MarkerFaceColor','k');
% %plot( sim_normed_time/1e6 + timeOffset, sim_normed_droplet_diameter,'Color','b','MarkerFaceColor','b');
% %plot( sim_normed_time/1e6 + timeOffset, sim_normed_droplet_diameter);
% xlabel('t/d_0^2 [s/mm^2]')
% ylabel('d^2/d_0^2 [-] ')
% xlim([0 1])
% legend('experiment','simulation');

% figure(2)
% %subplot(3,2,1)
% box on
% hold on
% plot( time(1:size_active)*1e6./(s_radius(1)*2*1e6).^2, sim_normed_droplet_diameter(1:size_active),'LineWidth',2,'Color','b' )
% %plot experiment
% xlabel('t/d_0^2 [s/mm^2]')
% ylabel('(d/d_0)^2 ')
% %legend('experiment','simulation');
% % xlim([0 4])
% % ylim([0 2])
% % 
% subplot(3,2,2)
% box on
% hold on
% plot( s_radius./s_radius(1), s_surfaceConcentration./s_surfaceConcentration(1));
% xlabel('R(t)/R_0')
% ylabel('particle concentration change N/N_0')
% set(gca, 'yscale', 'log')
% 
% subplot(3,2,3)
% box on
% hold on
% plot( time*1e3, s_volFrac);
% if lockingPoint > 0
%     plot( time(lockingPoint)*1e3, s_volFrac(lockingPoint),'Marker', 'o','MarkerFaceColor','r','Color','r');
% end
% xlabel('t / ms')
% ylabel('particle volume fraction \phi')
% set(gca, 'yscale', 'log')
% 
% subplot(3,2,4)
% box on
% hold on
% plot( time, s_massConservationCheck./s_massConservationCheck(1));
% xlabel('time / s')
% ylabel('total mass')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% subplot(3,2,5)
% box on
% hold on
% plot( time(1:length(s_particle_diam_surf))*1e3, s_particle_diam_surf*1e9);
% plot( time(1:length(s_particle_diam_surf))*1e3, s_particle_diam_core*1e9);
% set(gca,'yscale','log')
% xlabel('time / ms')
% ylabel('particle diameter / nm')
% 
% subplot(3,2,6)
% box on
% hold on
% %plot( time, s_collisionRate);
% %plot( time, s_collisionRate_diff);
% plot( time, s_collisionRate_reg);
% %plot( time, s_collisionRate_circ);
% set(gca,'yscale','log')
% xlabel('time / s')
% ylabel('collision rate / #/m³/s')
% %legend('collisionRate','collisionRate brownian','collisionRate regression','collisionRate circulation')
% 
% figure(3)
% box on;
% hold on;
% %plot(droplet_LifeTime, s_collisionRate_reg(1:size_active-1)./s_radius(1:size_active-1));
% plot(droplet_LifeTime, s_collisionRate_reg(1:size_active-1)./s_surfaceShellWidth(1:size_active-1));
% set(gca,'yscale','log')
% xlabel('1 - d^3/d_0^3')
% xlim([0 1])
% ylabel('\beta_{reg} / R ')


%% plot psds
meanDiam_core = meandiamg(gridDiams(3:end)*1e9, particleConc_core(3:end))
meanDiam_surf = meandiamg(gridDiams(3:end)*1e9, particleConc_surf(3:end));

particle_Conc_core_normed = normPSD(gridDiams*1e9, particleConc_core);
particle_Conc_surf_normed= normPSD(gridDiams*1e9, particleConc_surf);

volMeanDiam = volWeightedMeanDiam(gridDiams, particleConc_surf)*1e9;

figure(5)
subplot(2,1,1)
box on;
hold on;
%not normalized yet
plot(gridDiams(1:end-1)*1e9, particle_Conc_core_normed,'Color','b');
plot([meanDiam_core meanDiam_core], [1 max(particleConc_core)],'Color','b');
plot(gridDiams(1:end-1)*1e9, particle_Conc_surf_normed, 'Color','r' );
plot([meanDiam_surf meanDiam_surf], [1 max(particleConc_surf)],'Color','r');
plot([volMeanDiam volMeanDiam], [1 max(particleConc_surf)],'Color','r');
%plot vertical lines
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('particle diameter / nm')
ylabel('dN/dlog d_p / m³')
ylim([1 max(particleConc_surf + particleConc_core)])
legend('core','core mean diam','surface','surface mean diam')

particleVolumes_core = particleConc_core.*gridVols;
sumVolCore = sum(particleConc_core.*gridVols);
particleVolumes_surf = particleConc_surf.*gridVols;
sumVolSurf = sum(particleConc_surf.*gridVols);

particle_Vol_surf_normed= normPSD(gridDiams*1e9, particleVolumes_surf);

% figure(6)
% box on;
% hold on;
% %not normalized yet
% plot(gridDiams(1:end-1)*1e9, particle_Vol_surf_normed/sum(particle_Vol_surf_normed), 'Color','r' );
% %plot vertical lines
% set(gca,'xscale','log')
% set(gca,'yscale','lin')
% xlabel('particle diameter / nm')
% ylabel('dV/dlog d_p / V_{total} / m³')
% ylim([1e-6 1])

%psd in gas phase
V_surf = calcShellVolume(r, r - surfaceShellWidth);
V_core = calcShellVolume(r - surfaceShellWidth, 0);
sum(particleConc_core*V_core)*6e3
concentrationsGasPhase = particleConc_core*V_core + particleConc_surf*V_surf;
mean_diam_gas = meandiamg(gridDiams(2:end)*1e9, concentrationsGasPhase(2:end));
mean_diam_gas_core = meandiamg(gridDiams(2:end)*1e9, particleConc_core(2:end)*V_core);
noParticlesProduced = sum(concentrationsGasPhase(2:end));
noParticlesProduced_core = sum(particleConc_core(2:end)*V_core);
sigma_g = sigmag(gridDiams(2:end), concentrationsGasPhase(2:end));

concentrationsGasPhaseNormed = particle_Conc_core_normed*V_core + particle_Conc_surf_normed*V_surf;

% figure(7)
% title('Core particles')
% box on;
% hold on;
% plot(gridDiams(1:end-1)*1e9, particle_Conc_core_normed*V_core,'Color','b');
% plot([mean_diam_gas mean_diam_gas], [1 max( particle_Conc_core_normed(3:end)*V_core) ],'Color','b');
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('particle diameter / nm')
% ylabel('dN/dlog d_p / m³')
% ylim([1 max( particle_Conc_core_normed*V_core )])

mean_diam_gas_surf = meandiamg(gridDiams(2:end)*1e9, particleConc_surf(2:end)*V_surf);

% figure(8)
% title('Surface shell particles')
% box on;
% hold on;
% plot(gridDiams(1:end-1)*1e9, particle_Conc_surf_normed*V_surf,'Color','b');
% plot([mean_diam_gas_surf mean_diam_gas_surf], [1 max( particle_Conc_surf_normed(3:end)*V_surf) ],'Color','b');
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('particle diameter / nm')
% ylabel('dN/dlog d_p / m³')
% ylim([1 max( particle_Conc_surf_normed*V_surf )])

%with breakage model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particleConcAccumulated_normed = normPSD(gridDiams*1e9, particleConcAccumulated);
% mean_diam_breakage = meandiamg(gridDiams*1e9, particleConcAccumulated);
% figure(8)
% title('With Breakage model')
% box on;
% hold on;
% plot(gridDiams(1:end-1)*1e9, particleConcAccumulated_normed,'Color','b');
% plot([mean_diam_breakage mean_diam_breakage], [1 max( particleConcAccumulated_normed) ],'Color','b');
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('particle diameter / nm')
% ylabel('dN/dlog d_p / m³')
% ylim([1 max( particleConcAccumulated_normed )])

%plot comparison with Suleiman et al., 2021
% CMD_Suleiman = [4.45 5.63 7.01 7.81];
% CMD_calculated = [5.054 5.648 5.657 4.327];
% CMD_calculated_breakage = [5.054 5.648 6.273 6.133];
% conc = [0.05 0.1 0.2 0.5]
% 
% figure
% hold on
% box on
% plot(conc, CMD_Suleiman, 'LineWidth',2, 'Color','k','Marker','o','MarkerFaceColor','k', 'MarkerSize',4)
% plot(conc, CMD_calculated ,'LineWidth',2,'Color','b','Marker','d','MarkerFaceColor','b')
% plot(conc, CMD_calculated_breakage,'LineWidth',2,'Color','b','LineStyle',':','Marker','v','MarkerFaceColor','b')
% xlabel('precursor concentration [mol/l]')
% ylabel('d_{5,0} [nm]')
% ylim([0 10])
% legend('Suleimant et al, 2021','No Breakage','With Breakage')
% 
% figure
% hold on
% box on
% for i = 1:length(conc)
%     plot(conc(i), CMD_calculated(i))
% end
% xlabel('precursor concentration [mol/l]')
% ylabel('d_{5,0} [nm]')
% ylim([0 10])