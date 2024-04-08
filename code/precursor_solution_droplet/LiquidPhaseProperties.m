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
        calgToJkg = 4186.7999; % unit conversion factor
        calToJ = 4.186; % unit conversion factor
        R = 8.314; %J/mol/K, ideal gas constant 
    end
    
    methods
        function obj = LiquidPhaseProperties(compoundName)
            obj.name = compoundName;
                                           
            if strcmp(compoundName,'2EHA') 
            % Values taken from Kunstmann et al. (2023): 
            % "Thermophysical Properties of Mixtures of 
            % 2Ethylhexanoic Acid and Ethanol" for T = 333.15 K 
            % Here, we use pure 2EHA as an approximate for a 
            % mixture of 2EHA and ethanol, because EtOH quickly 
			% evaporates from the droplet surface during
            % combustion

                obj.M = 0.144214;
                obj.rho = 873.8; % kg/m³
                obj.T_boil = 501.2; % K
                obj.lambda_heat = 0.133; % W/m/K 
                obj.cp = 274/obj.M; %J/kg/K 
                obj.L = 4.2857e5; % J/kg at 298 K
                obj.Q = 3.3284e7; % J/kg
                obj.nu = 4.0969e-01;           
                obj.dyn_visc = 2.413e-3; % dynamic viscosity, Pa*s
            end
            
            obj.Rg = obj.R/obj.M;
            obj.D_heat = obj.lambda_heat/obj.rho/obj.cp;
        end

    end
end

