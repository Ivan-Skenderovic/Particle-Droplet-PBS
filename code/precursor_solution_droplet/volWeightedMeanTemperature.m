function Tw = volWeightedMeanTemperature(T, diameter)

noSpatialGridPoints = length(T);
dropletVolume = pi/6*diameter^3;
delta_r = diameter / 2 / noSpatialGridPoints;

dV_i = 0;
Tw = 0;
r_inner = 0;
r_outer = 0;

    for i = 0 : (noSpatialGridPoints - 1)
	
        r_outer = (i + 1)*delta_r;
        
        r_inner =  i*delta_r;
        
        dV_i = 4/3*pi*(r_outer^3 - r_inner^3)/dropletVolume;
        
        Tw = Tw + dV_i*T(i + 1);
		
    end

end

