function collRates = collisionRatesReg(coagConst, gridVolumes, noNodes)

collRates = zeros(noNodes, noNodes);
firstParticleNode = 3;

    for i= firstParticleNode : noNodes  
        for j = firstParticleNode : noNodes
            
            collRates(i, j) = coagConst* ...
          ( (3*gridVolumes(i)/4/pi)^(1/3) + (3*gridVolumes(j)/4/pi)^(1/3) )^2;            
        
        end
    end
    
end

