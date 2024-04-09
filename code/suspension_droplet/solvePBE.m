function dYdt = solvePBE(~ , N, gridVolumes, beta, splittingOps, reactionRate,...
    firstParticleNode)

noNodes = length(N);
dNdt(1:noNodes) = 0;

%% solve precursor decomposition numerically (example): %%%%%%%%%%%%%%%%%%%
% A reaction rate of type dCdt = -reactionRate*C is assumed, where C
% is the precursor concentration. In case the user wants to
% supply their own chemistry model, dNdt(1), dNdt(2) need to be set
% accordingly and criticalClusterSize and nucleationRate need to be defined
% according to the chosen nucleation theory. 

dNdt(1) = -N(1)*reactionRate; 
noMolecsInCriticalCluster = 1; 
nucleationRate = -dNdt(1)/noMolecsInCriticalCluster; 
dNdt(2) = -dNdt(1) - nucleationRate*noMolecsInCriticalCluster;

kelvinVolume = gridVolumes(2); % equal to monomer volume in this study
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
    
    dNdt(k) = nucleationRate*noMolecsInCriticalCluster*nucleationSplit(k) + ...
              0.5*inFlowSum - N(k)*outFlowSum; 
    
end

dYdt = dNdt';

end

