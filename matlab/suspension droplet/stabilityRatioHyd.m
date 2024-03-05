function W = stabilityRatioHyd(particleConc, gridVolumes, volFracCrit)
    % calculates the hydrodynamic interaction paramenter for a monodisperse
    % particle population according to Honig et al. (1971), "Effect of 
    % Hydrodynamic Interaction on the Coagulation Rate of Hydrophobic Colloids"
    % and calculated in the form given by Bal and Bandyopadhyaya (2018), 
    % "Generalized Model for Nano- and Submicron Particle Formation in
    % Liquid Phase, Incorporating Reaction Kinetics and Hydrodynamic Interaction: 
    % Experiment, Modeling, and Simulation" 

    volFrac = sum(particleConc(3:end).*gridVolumes);
    r_a = 2*(volFracCrit/volFrac)^(1/3); %dimensionless mean particle 
    % center-to-center distance
    
    if r_a > 2
        G_hyd = (6*(r_a - 2)^2 + 4*(r_a - 2))/(6*(r_a - 2)^2 + 13*(r_a - 2) + 2);     
    else
        G_hyd = 1; %touching condition
    end

    W = 1/G_hyd; %assuming no energy barrier between particles

end

