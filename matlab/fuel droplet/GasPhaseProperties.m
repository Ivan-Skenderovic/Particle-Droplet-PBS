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
           obj.T_inf = 1600; %K
           obj.p_inf = 1;
           obj.Y_0inf = 0.232;
           obj.M_air = 0.02897; %air, kg/mole
           obj.lambda_heat = 0.02638; %at 300K and 1 bar for air
           obj.C_g = 1465.379998; %J/kg/K, not normalized
           obj.Cg_Law = 0.35; %cal/gm/K, normalized according to C.K. Law
        end    
    end
end

