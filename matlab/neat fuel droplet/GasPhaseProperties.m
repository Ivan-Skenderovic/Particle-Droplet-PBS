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
    end
    
    methods
        function obj = GasPhaseProperties()
           obj.T_inf = 1600;
           obj.p_inf = 1;
           obj.Y_0inf = 0.232;
           obj.M_air = 0.02897; %air
           %obj.M_air = 0.028; %nitrogen
           obj.lambda_heat = 0.02638; %at 300K and 1 bar for air
           %obj.lambda_heat = 0.02597; %at 300K and 1 bar for nitrogen
           obj.C_g = 1465.379998; %J/kg/K
           obj.Cg_Law = 0.35; %cal/gm/K
        end    
    end
end

