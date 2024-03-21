function W = stabilityRatioHyd(particleConc, gridVolumes, volFracCrit)
    % calculates the hydrodynamic interaction paramenter for a monodisperse
    % particle population according to Honig et al. (1971), "Effect of 
    % Hydrodynamic Interaction on the Coagulation Rate of Hydrophobic Colloids"
    % and calculated in the form given by Bal and Bandyopadhyaya (2018), 
    % "Generalized Model for Nano- and Submicron Particle Formation in
    % Liquid Phase, Incorporating Reaction Kinetics and Hydrodynamic Interaction: 
    % Experiment, Modeling, and Simulation" 
    
    volFrac = sum(particleConc(3:end).*gridVolumes(3:end));

    %check if particles are present
    if (volFrac == 0)
        W = 1; 
        return;
    end
   
    r_a = 2*(volFracCrit/volFrac)^(1/3); % dimensionless average particle 
    % center-to-center distance
    
    if r_a > 2 % minimum distance is 2 times average particle radius
        G_hyd = (6*(r_a - 2)^2 + 4*(r_a - 2))/(6*(r_a - 2)^2 + 13*(r_a - 2) + 2);     
    else %touching condition
        G_hyd = 1; 
    end

    W = 1/G_hyd; % assuming no energy barrier between particles -> 
    % interaction integrals equals 1

end

