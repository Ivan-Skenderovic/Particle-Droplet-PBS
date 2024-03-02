function collRates = collisionRatesDiff(coagConst, gridVolumes, noNodes)

collRates = zeros(noNodes, noNodes);
firstParticleNode = 3;

    for i= firstParticleNode : noNodes  
        for j = firstParticleNode : noNodes
            
            collRates(i, j) = coagConst* ...
            ( 1/gridVolumes(i)^(1/3) + 1/gridVolumes(j)^(1/3) ) * ...
            ( gridVolumes(i)^(1/3) + gridVolumes(j)^(1/3) );            
        
        end
    end
    
end

