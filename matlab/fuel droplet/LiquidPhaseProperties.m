classdef LiquidPhaseProperties
    
    properties        
        name; %name string
        M; % Molar mass in kg/mol
        Rg; % specific gas constant R/M in J/kg/K
        D_heat; % heat diffusion coefficient in m²/s calculated from lambda_heat/rho/cp  
        rho; % density in kg/m³
        T_boil; % boiling Temperature in K
        lambda_heat; % thermal conductivity in W/m/K
        cp; % heat capacity in J/kg/K
        L; %latent heat of vaporization in J/kg
        Q; %heat of reaction per unit mass of fuel reacted in J/kg
        nu; %fuel oxidizer mass ratio for complete and stochiometric combustion    
        dyn_visc; %dynamic viscosity in Pa*s
        calgToJkg = 4186.7999;
        calToJ = 4.186;
        R = 8.314; %J/mol/K, ideal gas constant 
    end
    
    methods
        function obj = LiquidPhaseProperties(compoundName)
            obj.name = compoundName;
            
            if strcmp(compoundName,'Octane')           
                obj.M = 0.114; %kg/mol 
                obj.rho = 703; %kg/m³
                obj.T_boil = 399; %K
                obj.lambda_heat = 2.597e-01;%C.K. Law assumption W/m/K
                obj.cp = 1465.379998; %J/kg/K, C.K. Law  
                obj.L = 301449.59953; %J/kg, C.K. Law 
                obj.Q = 44798759.930132; %J/kg, C.K. Law 
                obj.nu = 0.284; %C.K. law
            end
            
            if strcmp(compoundName,'Ethanol')       
                obj.M = 0.04607;
                obj.rho = 789;
                obj.T_boil = 3.515200000000000e+02;
                obj.lambda_heat = 0.171;
                obj.cp = 112.4/obj.M;
                obj.L = 38560/obj.M;
                obj.Q = 1370700/obj.M;
                obj.nu = 0.479895;
                obj.dyn_visc = 0.4e-3; %Pa*s
            end
                     
            %average over all fuel species
            if strcmp(compoundName,'Mixture')
                obj.M = 0;
                obj.rho = 0;
                obj.T_boil = 0;
                obj.lambda_heat = 0;
                obj.cp = 0;
                obj.L = 0;
                obj.Q = 0;
                obj.nu = 0;
                obj.dyn_visc = 0.0; %Pa*s
            end
            
            obj.Rg = obj.R/obj.M;
            obj.D_heat = obj.lambda_heat/obj.rho/obj.cp;
        end

    end
end

