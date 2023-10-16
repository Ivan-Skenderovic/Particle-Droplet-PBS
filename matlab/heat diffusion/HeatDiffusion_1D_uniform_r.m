function [Y] = HeatDiffusion_1D_uniform_r(t, TandDiam, boundaryType, LiquidData, GasData)

ConfigSettings

noTemperatures = length(TandDiam) - 1;
Diam2 = TandDiam(end);

if(sqrt(Diam2) > 1e-7)    
   % r = sqrt(Diam2)/2;
    r = 1;
    LiquidData.D_heat = LiquidData.D_heat/(Diam2/4);
    delta_r = r/NOSPATIALGRIDPOINTS;
    if ~isreal(delta_r)
        warning('delta_r');
    end
        
    r(1:noTemperatures) = delta_r:delta_r:r;

    dTdt(1:noTemperatures) = 0;
    innerNodeIdx = 2:noTemperatures - 1;

    % center temperature node
    dTdt(1) = LiquidData.D_heat*6/delta_r^2*(-TandDiam(1) + TandDiam(2));

    %inner temperature nodes
    dTdt(innerNodeIdx) = LiquidData.D_heat*( ...
    (TandDiam(innerNodeIdx-1) - 2.*TandDiam(innerNodeIdx) + TandDiam(innerNodeIdx+1))/delta_r.^2 ...
    + 1./r(innerNodeIdx)'*1./delta_r.*(TandDiam(innerNodeIdx+1) - TandDiam(innerNodeIdx-1)));

    %surface temperature node
    if boundaryType == 1
        %insulation    
    elseif boundaryType == 2
        %temperature dependent heat flux
            dTdt(noTemperatures) = 2*LiquidData.D_heat/delta_r^2*(TandDiam(noTemperatures-1)...
                - TandDiam(noTemperatures)...
                -(1+1/noTemperatures)*LiquidData.Q*delta_r*(TandDiam(noTemperatures) - GasData.T_inf));
    elseif boundaryType == 3
        %sirignano boundary conditions for droplet evaporation
        Rg = LiquidData.Rg/LiquidData.calgToJkg; %cal/gm
        %specific latent heat of vaporization
        L = LiquidData.L/LiquidData.calgToJkg; %cal/gm
        %thermal conductivity of gas phase
        lambda_heat_g = GasData.lambda_heat;
        %thermal conductivity of liquid phase
        lambda_heat_l = LiquidData.lambda_heat; %W/m/K
        %specific gas and liquids latent heats 
        c_p = LiquidData.cp/LiquidData.calgToJkg; %cal/gm
        %heat of reaction per unit mass of fuel reacted
        Q = LiquidData.Q/LiquidData.calgToJkg; %cal/gm
        %stochiometric fuel oxidizer mass ratio
        nu = LiquidData.nu;
        %boiling point
        T_boil = LiquidData.T_boil; %K
        p_inf = GasData.p_inf; %atm
        %oxidizer mass fraction
        Y_0inf = GasData.Y_0inf;
        %ambient temperature
        Tinf = GasData.T_inf; 
        T_dinf = c_p/L*GasData.T_inf; 
        T_ds = c_p/L.*TandDiam(noTemperatures);
        T_db = c_p/L*T_boil;

        M_air = GasData.M_air; %kg/mol
        W_A = M_air; %average of all species except fuel in kg/mol
        W_F = LiquidData.M; %kg/mol

        Y_Fs = (1 + (W_A/W_F)*(p_inf*exp((c_p/Rg)*(1./T_ds - 1/T_db))-1)).^-1;
        if Y_Fs > 1
            warning('Y_Fs');
        end
        %dimensionless surface temperature        
        H_d = (1 - Y_Fs)*(T_dinf - T_ds + nu*Y_0inf*Q/L) / (Y_Fs + nu*Y_0inf);
        
       % B = (nu*Y_0inf + Y_Fs)/(1 - Y_Fs);
        B = (T_dinf - T_ds*nu*Y_0inf*Q/L)/H_d;
        
        if B < -1
            warning('B');
        end
        
        A = log(1+B)./r(noTemperatures)*lambda_heat_g./lambda_heat_l;
        A_= A/B*Tinf + A/B*nu*Q*Y_0inf/c_p - A*L/c_p;

        dTdt(noTemperatures) = LiquidData.D_heat/delta_r^2*(1-1/noTemperatures)*TandDiam(noTemperatures-1)...
            -2*LiquidData.D_heat/delta_r^2*TandDiam(noTemperatures)...
            + LiquidData.D_heat/delta_r^2*(1+1/noTemperatures)*( 2*delta_r*(A_ - A/B*TandDiam(noTemperatures))...
            +  TandDiam(noTemperatures-1));
    end
    
    %dimensionless evoration rate
    m_d = log(1 + (T_dinf - T_ds + nu*Y_0inf*Q/L)/H_d);
    %correction for droplet evaporation in convective flow
    Re_Droplet = Air_DENSITY*sqrt(Diam2) * DROPLET_RELATIVE_VELOCITY / Air_DYNAMIC_VISCOSITY;
    m_d = m_d*(1 + 0.3*Air_PRANDTL^(1.0/3.0)*(2*Re_Droplet)^0.5);
    %dr²dt = -K*t calculated here, dd²dt = 4*dr²dt needed for beta0
    surfRegRate = 4*m_d * 2 * lambda_heat_g / Ethanol_DENSITY / Air_HEAT_CAP;

    dDiam2dt = - surfRegRate;
    dTdt = dTdt';

    Y = [dTdt; dDiam2dt];
else

dTdt(1:noTemperatures) = 0;
dDiam2dt = 0;
   
Y = [dTdt'; dDiam2dt];    

end

