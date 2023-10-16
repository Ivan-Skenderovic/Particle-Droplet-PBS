classdef LiquidPhaseProperties
    
    properties        
        name; %name string
        C; %No of carbon atoms in molecule
        H; %No of hydrogen atoms in molecule
        M; % Molar mass in kg/mol
        M_anhydrous; % anhydrous Molar mass
        Rg; % specific gas constant R/M in J/kg/K
        D_heat; % heat diffusion coefficient in m²/s calculated from lambda_heat/rho/cp  
        D_g; %mass diffusion coefficient in Air/N2 m²/s
        D_self; %Self Diffusion coefficient m²/s
        rho; % density in kg/m³
        T_boil; % boiling Temperature in K
        p_boil; %boiling pressure in Pa
        lambda_heat; % thermal conductivity in W/m/K
        cp; % heat capacity in J/kg/K
        L; %latent heat of vaporization in J/kg
        Q; %heat of reaction per unit mass of fuel reacted in J/kg
        nu; %fuel oxidizer mass ratio for complete and stochiometric combustion
            %calculated from combustion equation O2
        X_liq; %Molar fraction of species at droplet surface liquid side, default value is 1
        X_gas; %Molar fraction of species at droplet surface gas side, default value is 1
        Y_liq; %Mass fraction of species in fuel, default value is 1
        Y_gas; %Mass fraction of species in fuel, default value is 1
        epsilon; %fractional evaporation rate m_i/m according to Law model
        Y_h2o_in_Prec;
        %Antoine Parameter:
        A_Antoine;
        B_Antoine;
        C_Antoine;        
        sol_aq; %solubility in water, mol/l
        sol_ethanol; %solubility in ethanol, mol/l
        dyn_visc; %dynamic viscosity in Pa*s
        r_hyd; %hydrodynamic radius of molecule
        initial_diameter; %initial diameter of droplet
        D_Stokes; %diffusion coefficient in liquids according to einstein-
        %stokes relation in m²/s
        p_crit; %critical pressure
        T_crit; %critical temperature
        T_superheat; %temperature limit at superheat occurs
        calgToJkg = 4186.7999;%converts values between cal/g and J/kg
        calToJ = 4.186; %converts values between cal and J
        R = 8.314; %ideal gas constant J/mol/K
        NA = 6.022*10^23; %1/mol
        kB = 1.30*10^-23;
        sol_anhydrous; % 0 or 1, depending if soluble in anhydrous form
        plotColor;
    end
    
    methods
        function obj = LiquidPhaseProperties(compoundName)
            
            obj.name = compoundName;
            
            %C8_H_18
            if strcmp(compoundName,'Octane')           
                obj.M = 0.114; %kg/mol 
                obj.rho = 703; %kg/m³
                obj.T_boil = 399; %K
                %obj.lambda_heat = 0.137; % W/m/K
                obj.lambda_heat = 2.597e-01;%C.K. law assumption W/m/K
                %obj.cp = 255.68/obj.M; %J/kg
                obj.cp = 1465.379998; %C.K. law assumption J/kg/K
                obj.L = 301449.59953; %C.K. law J/kg
                obj.Q = 44798759.930132; %C.K. law J/kg
                obj.nu = 0.284; %C.K. law
            end
            
            %C2_H6_OH
            if strcmp(compoundName,'Ethanol')    
            obj.name = 'Ethanol';       
            obj.M = 0.04607;
            obj.rho = 789;
            obj.T_boil = 3.515200000000000e+02;
            obj.lambda_heat = 0.171;
            obj.cp = 112.4/obj.M;
            obj.L = 38560/obj.M;
            obj.Q = 1370700/obj.M;
            obj.nu = 0.479895;
            obj.p_crit = 63.06439; %critical pressure in atm
            obj.T_crit = 516.25; %critical temperature
            obj.dyn_visc = 0.4e-3; %Pas
            %obj.T_superheat = 0; %temperature limit at superheat occurs
            end
            
            %C4_H10_O
            if strcmp(compoundName,'Butanol')  
            obj.name = 'Butanol';
            obj.M = 0.074123;
            obj.rho = 809.5;
            obj.T_boil = 390.88;
            obj.lambda_heat = 0.1499;
            obj.cp = 192.2/obj.M;
            obj.L = 43.29e3/obj.M;
            obj.Q = 2676600/obj.M;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0; 
            obj.epsilon = 0;
            obj.nu = 0.386057;
            end
            
            %C7_H16
            if strcmp(compoundName,'Heptane')
            obj.name = 'Heptane';
            obj.M = 0.10021;
            obj.D_g = 6.54404991*1e-6;
            obj.D_self = 3.22e-9;
            obj.rho = 679.5;
            obj.T_boil = 371.53;
            obj.lambda_heat = 0.15;
            obj.cp = 224.64/obj.M;
            obj.L = 31770/obj.M;
            obj.Q = 4.85e6/obj.M;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0; 
            obj.epsilon = 0; 
            obj.C = 7;
            obj.H = 16;
            obj.nu = 0;
            end
            
            %C16H34
            if strcmp(compoundName,'Hexadecane')
            obj.name = 'Hexadecane';
            obj.M = 0.22645;
            obj.D_g = 6.54404991*1e-6;
            obj.D_self = 0.378e-9;
            obj.rho = 770;
            obj.T_boil = 560;
            obj.lambda_heat = 0.15;
            obj.cp = 499.72/obj.M;
            obj.L = 59800/obj.M;
            obj.Q = 10.7e6/obj.M;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0; 
            obj.epsilon = 0; 
            obj.C = 0;
            obj.H =0;
            obj.nu = 0;
            end
            
            %C10H22
            if strcmp(compoundName,'Decane')
            obj.name = 'Decane';
            obj.M = 0.14229;
            obj.D_g = 5.42e-6;
            obj.D_self = 1.56e-9;
            obj.rho = 730;
            obj.T_boil = 447.15;
            obj.lambda_heat = 0.1168;
            obj.cp = 315.45/obj.M;
            obj.L = 51500/obj.M;
            obj.Q = 6779210/obj.M;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0; 
            obj.epsilon = 0; 
            obj.C = 10;
            obj.H = 22;
            %Antoine Parameter from 367-448 K, pressure in bar
            obj.A_Antoine = 4.07857;
            obj.B_Antoine = 1501.268;
            obj.C_Antoine = -78.67;  
            obj.nu = 0.269488;
            end
            
            %C12H26
            if strcmp(compoundName,'Dodecane')   
            obj.name = 'Dodecane';
            obj.M = 0.17034;
            obj.D_g = 4.81e-6;
            obj.D_self = 1.34e-9;
            obj.rho = 749.5;
            obj.T_boil = 491;
            obj.lambda_heat = 0.1123;
            obj.cp = 375.6/obj.M;
            obj.L = 51600/obj.M;
            obj.Q = 7901740/obj.M;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0; 
            obj.epsilon = 0;
            obj.C = 12;
            obj.H = 26;
            %Antoine Parameter from 400-490 K, pressure in bar
            obj.A_Antoine = 4.10549;
            obj.B_Antoine = 1625.928;
            obj.C_Antoine = -92.839;  
            obj.nu = 0.39430555;
            end
            
            %H2O
            if strcmp(compoundName,'Water') 
            obj.name = 'Water';
            obj.M = 0.018015268;
            obj.D_g = 0.282e-4;
            obj.D_self = 8.667e-9;
            obj.rho = 998;
            obj.T_boil = 372.34;
            %obj.T_boil = 500;
            obj.lambda_heat = 0.591;
            obj.cp = 4187; 
            obj.L = 40600/obj.M;
            obj.dyn_visc = 0.0005474;
            obj.Q = 0;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0;
            obj.epsilon = 0;
            obj.C = 0;
            obj.H = 2;
            obj.nu = 0;
            obj.p_crit = 217.75; %critical pressure
            obj.T_crit = 647; %critical temperature
            %obj.T_superheat= 0; %temperature limit at superheat occurs
            end
            
            %Al(NO3)3
            if strcmp(compoundName,'AluminiumIIINitrate') 
            obj.name = 'AluminiumIIINitrate';
            obj.M_anhydrous = 0.21299;
            obj.M = 0.37513; %as nonahydrate
            obj.D_g = 0;
            obj.D_self = 0;
            obj.rho = 1720;
            obj.T_boil = 423;
            obj.lambda_heat = 0.661;
            obj.cp = 3194; 
            obj.L = 0;
            obj.Q = 0;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0;
            obj.epsilon = 0;
            obj.C = 0;
            obj.H = 18;
            obj.nu = 0;
            obj.sol_aq = 0.673; %kg_nitrate nonahydrate /kg_water at 25°C,
                                %1.6 as  kg_nitrate anhydrous/kg_water
            obj.sol_ethanol = 0.107769; %kg_nitrate/kg_ethanol at 25°C
            obj.sol_anhydrous = 1;
            obj.p_crit = 0; %critical pressure
            obj.T_crit = 0; %critical temperature
            obj.dyn_visc = 0.891e-3;
            obj.Y_h2o_in_Prec = 0.43; % 43% water per gramm precursor
            %obj.T_superheat= 0; %temperature limit at superheat occurs
            obj.plotColor = [0 1 0];
            obj.initial_diameter = 90e-6;
            end
            
            %Fe(NO3)3
            if strcmp(compoundName,'IronIIINitrate') 
            obj.name = 'IronIIINitrate';
            obj.M = 0.40399;% as nonahydrate
            obj.M_anhydrous = 0.24186;
            obj.D_g = 0;
            obj.D_self = 0;
            obj.rho = 1680;
            obj.T_boil = 398;
            obj.lambda_heat = 0.594;
            obj.cp = 3215.1; 
            obj.L = 0;
            obj.Q = 0;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0;
            obj.epsilon = 0;
            obj.C = 0;
            obj.H = 0;
            obj.nu = 0;
            obj.sol_aq = 1.5; %kg_nitrate/kg_water, hexahydrate
            obj.sol_ethanol = 1/21; %kg_nitrate/kg_ethanol at 25°C
            obj.p_crit = 0; %critical pressure
            obj.T_crit = 0; %critical temperature
            obj.dyn_visc = 0.891e-3;
            obj.Y_h2o_in_Prec = 0.4; %  water per gramm precursor
            %obj.T_superheat= 0; %temperature limit at superheat occurs
            obj.plotColor = [1 0 0];
            obj.initial_diameter = 95e-6;
            end
            
            %Zn(NO3)2
            if strcmp(compoundName,'ZincIINitrate') 
            obj.name = 'ZincIINitrate';
            obj.M = 0.29749;% as hexahydrate
            obj.M_anhydrous = 0.18936; 
            obj.D_g = 0;
            obj.D_self = 0;
            obj.rho = 2065; %hexahydrate
            obj.T_boil = 398;
            obj.lambda_heat = 0.661;
            obj.cp = 2666.2; 
            obj.L = 0;
            obj.Q = 0;
            obj.X_liq = 0; 
            obj.Y_liq = 0; 
            obj.X_gas = 0; 
            obj.Y_gas = 0;
            obj.epsilon = 0;
            obj.C = 0;
            obj.H = 0;
            obj.nu = 0;
            obj.sol_aq = 1.843; %kg_nitrate/kg_water, as hexahydrate %3.27 as trihydrate
            obj.sol_ethanol = 1; %kg_nitrate/kg_ethanol at 25°C
            obj.p_crit = 0; %critical pressure
            obj.T_crit = 0; %critical temperature
            obj.dyn_visc = 0.891e-3;
            obj.Y_h2o_in_Prec = 0.36; % water percent per gramm precursor
            %obj.T_superheat= 0; %temperature limit at superheat occurs
            obj.plotColor = [0 0 1];
            obj.initial_diameter = 139e-6;
            end
            
            %average over all fuel species
            if strcmp(compoundName,'Average')
            %this values are initialized in the main program
            obj.M = 0;
            obj.D_g = 0;
            obj.D_self = 0;
            obj.rho = 0;
            obj.T_boil = 0;
            obj.lambda_heat = 0;
            obj.cp = 0;
            obj.L = 0;
            obj.Q = 0;
            obj.nu = 0;
            obj.p_crit = 0; %critical pressure
            obj.T_crit = 0; %critical temperature
            obj.dyn_visc = 0.0;
            obj.r_hyd = 0;
            %obj.T_superheat= 0; %temperature limit at 
            end
            
            obj.Rg = obj.R/obj.M;
            obj.D_heat = obj.lambda_heat/obj.rho/obj.cp;
            obj.r_hyd = 0.5*(obj.M/obj.NA/obj.rho*6/pi)^(1/3);
%           obj.D_Stokes = obj.kB*obj.T_boil/(6*pi*obj.dyn_visc*obj.r_hyd); 
        end
        
        function vapPressure = vapPressureAntoine(obj, Temperature)
            %in bar
            vapPressure = 10^(obj.A_Antoine - obj.B_Antoine/(Temperature + obj.C_Antoine)); 
        end
    end
end

