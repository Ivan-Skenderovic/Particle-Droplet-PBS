function W = stabilityRatioHyd(particleConc, gridDiams, volFracCrit)
    
    if sum(particleConc(3:end)) > 0
        particle_diam = meandiamg(gridDiams(3:end), particleConc(3:end));  
    else
        particle_diam = (6*particle_volume/pi)^(1/3);
    end
        
    noParticlesShell = sum(particleConc(3:end));
    
    volFrac = sum(particleConc(3:end).*pi/6.*gridDiams(3:end).^3);
    
    r_a = 2*(volFracCrit/volFrac)^(1/3);
    %r_a = 1/particle_diam*(3/(4*pi*noParticlesShell))^(1/3);
    %r_a = 1/particle_diam/noParticlesShell^(1/3);
    
    if r_a > 2
        G_hyd = (6*(r_a - 2)^2 + 4*(r_a - 2))/( 6*(r_a - 2)^2 + 13*(r_a - 2) + 2);     
    else
        G_hyd = 1; %touching condition
    end

    W = 1/G_hyd;

end

