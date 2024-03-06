function dYdt = nodal_ode(t, N, gridVolumes, beta, splittingOps, reactionRate)
noNodes = length(N);

criticalClusterSize = 1; % metal oxide negligible vapor pressure

dNdt(1:noNodes) = 0;

%precursor thermal decomposition
dNdt(1) = -N(1)*reactionRate;
% influx by decomposition and outflow by nucleation:
nucleationRate = -dNdt(1)/criticalClusterSize; %reaction rate controlled nucleation:
dNdt(2) = -dNdt(1) - nucleationRate*criticalClusterSize; 

firstParticleNode = 3;
kelvinVolume = gridVolumes(2); %equal to monomer
nucleationSplit(1:noNodes) = 0;

for k = firstParticleNode : noNodes
    
    outFlowSum = 0;
    inFlowSum = 0;

    %calculate nucleation split
    if kelvinVolume < gridVolumes(firstParticleNode)
        nucleationSplit(firstParticleNode) = ...
            kelvinVolume/gridVolumes(firstParticleNode);
    elseif gridVolumes(k-1) <= kelvinVolume && kelvinVolume < gridVolumes(k)
        nucleationSplit(k) = kelvinVolume/gridVolumes(k);           
    end
    
    for i = firstParticleNode : noNodes 
        outFlowSum = outFlowSum + beta(k,i)*N(i);
        for j = firstParticleNode : k      
           inFlowSum = inFlowSum + splittingOps(i,j,k).*beta(i,j).*...
                N(i)*N(j);            
        end        
    end
    
    dNdt(k) = nucleationRate*criticalClusterSize*nucleationSplit(k) + ...
              0.5*inFlowSum - N(k)*outFlowSum; 
    
end

dYdt = dNdt';

end

