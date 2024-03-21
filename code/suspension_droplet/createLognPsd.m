function particleConcentrations = createLognPsd(NONODES, gridDiams, PARTICLE_INITIAL_DIAMETER, ...
    sigma_g, ~, firstParticleNode, discretizationError)

err = 1;
particleConcentrations(1:NONODES) = 0;
particleMeanDiam = PARTICLE_INITIAL_DIAMETER;
counter = 0;

    while err > discretizationError
        for nodeIdx = firstParticleNode : NONODES
            
          particleConcentrations(nodeIdx) =  1/(gridDiams(nodeIdx)*log(sigma_g)*sqrt(2*pi)).* ...
                exp(-( log(particleMeanDiam) - log(gridDiams(nodeIdx)) ).^2 ...
                /2/log(sigma_g)^2 );
        end

        particleGeoMean_new = lognMedianDiam(gridDiams, particleConcentrations);

        err = abs((PARTICLE_INITIAL_DIAMETER - particleGeoMean_new)/PARTICLE_INITIAL_DIAMETER);

        particleMeanDiam = particleMeanDiam*(1 + discretizationError);

        counter = counter + 1;

        if counter > 1/discretizationError
            display('failed to meet tolerance in createLognPSD!');
            break
        end
    end
    
end
