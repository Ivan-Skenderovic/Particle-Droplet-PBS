function particleConcentrations = createLognPsd(NONODES, gridDiams, PARTICLE_INITIAL_DIAMETER, ...
    sigma_g, ~, firstParticleNode, discretizationError)
%Due to the particle volume discretization of the pivot method,
%an inherent error is introduced. To set the predefined moments of the 
%lognormal distribution with the mean diameter is shifted until the set
%moments are matched with the accuracy set by discretizationError.

err = 1;
particleConcentrations(1:NONODES) = 0;
%arbitrary starting value
particleMeanDiamShifted = PARTICLE_INITIAL_DIAMETER*0.5; 
counter = 0;

    while err > discretizationError
        for nodeIdx = firstParticleNode : NONODES            
		  %a lognormal distribution is assumed
          particleConcentrations(nodeIdx) =  1/(gridDiams(nodeIdx)*log(sigma_g)*sqrt(2*pi)).* ...
                exp(-( log(gridDiams(nodeIdx)) - log(particleMeanDiamShifted) ).^2 ...
                /2/log(sigma_g)^2 );
        end

        particleGeoMean_new = weightedGeoMean(gridDiams, particleConcentrations);
        err = abs((PARTICLE_INITIAL_DIAMETER - particleGeoMean_new)/PARTICLE_INITIAL_DIAMETER);

        particleMeanDiamShifted = particleMeanDiamShifted*(1 + discretizationError);

        counter = counter + 1;

        if counter > 1/discretizationError*1e3
            display('failed to meet tolerance in createLognPSD!');
            break
        end
    end
    
end
