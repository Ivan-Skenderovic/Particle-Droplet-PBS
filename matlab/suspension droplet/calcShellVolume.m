function V_shell = calcShellVolume(R_outer, r_inner)
    
  V_shell = 4./3*pi*(R_outer.^3 - r_inner.^3);    

end

