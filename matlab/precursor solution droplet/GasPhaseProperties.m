classdef GasPhaseProperties
   
    properties
        T_inf; %gas temperature in K
        p_inf; %ambient pressure in atm
        Y_0inf; % oxidizer mass fraction in gas phase
        M_air;  %kg/mol
        lambda_heat;  %heat conductivity in W/m/K
        rho_N2; %kg/m³
        C_g; %J/kg/K
        Cg_Law; %cal/gm
        X_air; %molar fraction of oxidizer gas at stochimetric combustion
        Y_air; %mass fraction of oxidizer gas at stochimetric combustion
        D_g; %self diffusion coefficient in m²/s
        calgToJkg = 4186.7999; % unit conversion factor
    end
    
    methods
        function obj = GasPhaseProperties()
           obj.T_inf = 2300; % K
           obj.p_inf = 1;
           obj.Y_0inf = 0.232;
           obj.M_air = 0.02897; % air, kg/mole
           obj.lambda_heat = 0.1175; % W/m/K, air at 2000 K and 0.1 MPa 
           obj.C_g = 1249.9; % J/kg/K, air at 1873 K and 0.1 MPa 
           obj.Cg_Law = 1249.9/obj.calgToJkg; % in cal/gm/K, as used by C.K. Law
        end    
    end
end

