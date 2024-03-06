function volumes = nodalGridVolumes(precursorMonomerVolume, monomerVolume, ...
    smallestParticleVolume, noNodes, gridSpacingFactor, firstParticleNode)

    qSpacingFactor = 10.0^(gridSpacingFactor / (noNodes - 2));

	volumes(1) = precursorMonomerVolume;
    volumes(2) = monomerVolume;
    
    if noNodes > 2
        for i = firstParticleNode : noNodes
            
            %original: qSpacingFactor^(i-1)
            volumes(i) = smallestParticleVolume * qSpacingFactor^(i - 2);

        end    
    end
    
end

