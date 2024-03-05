function collRates = collisionRatesDiff_fmr(coagConst, gridVolumes, noNodes)
% collision frequency function for the free molecular regime

collRates = zeros(noNodes, noNodes);
firstParticleNode = 3;

    for i= firstParticleNode : noNodes  
        for j = firstParticleNode : noNodes

             collRates(i, j) = coagConst* ...
               sqrt(1./gridVolumes(i) + 1./gridVolumes(j))*...
               (gridVolumes(i)^(1/3) + gridVolumes(j)^(1/3)).^2 ;        
        
        end
    end
    
end

