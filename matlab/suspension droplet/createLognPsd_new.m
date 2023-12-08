function [particleConcentrations, gridVols_new, indexPivotPoint] = ...
    createLognPsd(NONODES, gridVols, PARTICLE_INITIAL_DIAMETER, ...
    PARTICLE_INITIAL_CONCENTRATION, sigma_g, discretizationError)

gridDiams = (6.*gridVols./pi).^(1/3);
particleConcentrations(1:NONODES) = 0;
firstParticleNode = 3;

%set initial values for iteration
err = 1; 
particleGeoMean = PARTICLE_INITIAL_DIAMETER;
counter = 0;

    if sigma_g > 1
    % initialize a lognormal size distribution
        while err > discretizationError
            for nodeIdx = firstParticleNode : NONODES
              particleConcentrations(nodeIdx) =  ...
                    1/(gridDiams(nodeIdx)*log(sigma_g)*sqrt(2*pi)).* ...
                    exp(-( log(particleGeoMean) - log(gridDiams(nodeIdx)) ).^2 ...
                    /2/log(sigma_g)^2 );
            end
        
            particleGeoMean_new = meandiamg(gridDiams, particleConcentrations);        
            err = abs((PARTICLE_INITIAL_DIAMETER - particleGeoMean_new)/PARTICLE_INITIAL_DIAMETER);
            particleGeoMean = particleGeoMean*(1 + discretizationError);        
            counter = counter + 1;
        
            if counter > 1/discretizationError
                display('failed to meet tolerance in createLognPSD!');
                break
            end
        end
        
        indexPivotPoint = [];

        %rescale psd to set PARTICLE_INITIAL_CONCENTRATION
        particleConcentrations = particleConcentrations.*...
            PARTICLE_INITIAL_CONCENTRATION/sum(particleConcentrations);

    else
    % initialize a monodisperse population by replacing a pivot point with
    % size PARTICLE_INITIAL_DIAMETER
        indexPivotPoint = find(PARTICLE_INITIAL_DIAMETER < gridDiams,1);        
        gridVols(indexPivotPoint) = pi/6.*PARTICLE_INITIAL_DIAMETER^3;                
        particleConcentrations(indexPivotPoint) = PARTICLE_INITIAL_CONCENTRATION;
    
    end

    gridVols_new = gridVols;
end
