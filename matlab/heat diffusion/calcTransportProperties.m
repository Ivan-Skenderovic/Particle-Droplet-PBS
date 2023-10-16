%transport properties 
       
    ConfigSettings
    %sirignano boundary conditions for droplet evaporation
    Rg = LiqData.Rg/LiqData.calgToJkg; %cal/gm
    %specific latent heat of vaporization
    L = LiqData.L/LiqData.calgToJkg; %cal/gm
    %thermal conductivity of gas phase
    lambda_heat_g = GasData.lambda_heat;
    %thermal conductivity of liquid phase    
    lambda_heat_l = LiqData.lambda_heat;
    %specific gas and liquids latent heats 
    %c_p = LiqData.cp/LiqData.calgToJkg; %cal/gm
    c_p = GasData.Cg_Law; %cal/gm
    c_p_Air = GasData.C_g; % J/kg/K
    %heat of reaction per unit mass of fuel reacted
    Q = LiqData.Q/LiqData.calgToJkg; %cal/gm
    %stochiometric fuel oxidizer mass ratio
    nu = LiqData.nu;
    %boiling point
    T_boil = LiqData.T_boil; %K
    p_inf = GasData.p_inf; %atm
    %oxidizer mass fraction
    Y_0inf = GasData.Y_0inf;
    YFf = 0; % complete combustion
    %ambient temperature
    T_dinf = c_p/L*GasData.T_inf; 
    T_ds = c_p/L.*surfaceTemperature;
    T_db = c_p/L.*T_boil;

    M_air = GasData.M_air; %kg/mol
    W_A = M_air; %average of all species except fuel in kg/mol
    W_F = LiqData.M; %kg/mol

    Y_Fs = (1 + (W_A/W_F)*(p_inf*exp((c_p/Rg)*(1./T_ds - 1./T_db))-1)).^-1;  
    
    H_d = (1 - Y_Fs).*(T_dinf - T_ds + nu.*Y_0inf.*Q./L) ./ (Y_Fs + nu*Y_0inf - YFf.*(1 + nu*Y_0inf));
    
   %dimensionless evoration rate
    m_d = log(1 + (T_dinf - T_ds + nu.*Y_0inf.*Q./L)./H_d);    
    
    %correction for droplet evaporation in convective flow
    %Re_Droplet = Air_DENSITY*2*r(end)* DROPLET_RELATIVE_VELOCITY / Air_DYNAMIC_VISCOSITY;
    %m_d = m_d*(1 + 0.3.*Air_PRANDTL.^(1.0/3.0)*(2.*Re_Droplet).^0.5);
    