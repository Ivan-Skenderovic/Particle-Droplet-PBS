s_surfaceConcentration(i) = sum(particleConc_surf_new);
s_volFrac(i) = volFrac_surf;
s_surfaceShellWidth(i) = surfaceShellWidth;
s_radius(i) = r;   
s_particle_diam_surf_scale(i) = particle_diam_surf_scale;
s_particle_diam_core(i) = particle_diam_core;
s_meandiamg(i) = meandiamg(gridDiams(3:end), particleConc_core(3:end));
s_massConservationCheck(i) = massConservationError;
s_massConservationMergeCheck(i) = massConsErrorMerge;
s_massConservationAdvectionCheck(i) = massConsErrorAdvection;
s_massConservationSolverCheck(i) = massConsErrorSolver; 
s_concLastNode_surf(i) = particleConc_surf(end);
s_concLastNode_core(i) = particleConc_core(end);
   
