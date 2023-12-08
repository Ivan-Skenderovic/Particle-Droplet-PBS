function [Y] = HeatDiffusion_1D(t, TandRsquare, LiqData, GasData)

Rsquare = TandRsquare(end);
r = sqrt(Rsquare);
noTemperatures = length(TandRsquare) - 1;

if r < 1e-9 % terminate simulation at droplet size limit
    dTdt(1:noTemperatures) = 0;
    dr2dt = 0;
    Y = [dTdt'; dr2dt];   
    return
end

Temperatures = TandRsquare(1:end-1);
delta_r = r/noTemperatures;    
r(1:noTemperatures+1) = 0: delta_r: r;
dTdt(1:noTemperatures) = 0;
innerNodeIdx = 2:noTemperatures - 1;

%center temperature node
dTdt(1) = LiqData.D_heat*6/delta_r^2*(-Temperatures(1) + Temperatures(2));

%inner temperature nodes
dTdt(innerNodeIdx) = LiqData.D_heat*( ...
(Temperatures(innerNodeIdx-1) - 2.*Temperatures(innerNodeIdx) + Temperatures(innerNodeIdx+1))/delta_r.^2 ...
+ 1./r(innerNodeIdx)'*1./delta_r.*(Temperatures(innerNodeIdx+1) - Temperatures(innerNodeIdx-1)));

surfaceTemperature = Temperatures(end);
    
ConfigSettings

Rg = LiqData.Rg/LiqData.calgToJkg; %cal/gm
%specific latent heat of vaporization
L = LiqData.L/LiqData.calgToJkg; %cal/gm
%thermal conductivity of gas phase
lambda_heat_g = GasData.lambda_heat;
%thermal conductivity of liquid phase    
lambda_heat_l = LiqData.lambda_heat;
%specific gas and liquids latent heats 
c_p = GasData.Cg_Law; %cal/gm
c_p_Air = GasData.C_g; %J/kg/K
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

%Fuel mass fraction at droplet surface on gas side
Y_Fs = (1 + (W_A/W_F)*(p_inf*exp((c_p/Rg)*(1./T_ds - 1./T_db))-1)).^-1; 

H_d = (1 - Y_Fs).*(T_dinf - T_ds + nu.*Y_0inf.*Q./L) ./ ...
    (Y_Fs + nu*Y_0inf - YFf.*(1 + nu*Y_0inf));

%dimensionless evoration rate
m_d = log(1 + (T_dinf - T_ds + nu.*Y_0inf.*Q./L)./H_d);    

%Spalding heat transfer number
B = (T_dinf - T_ds + nu*Y_0inf*Q/L)/H_d;        
A = log(1+B)./r(noTemperatures)*lambda_heat_g/lambda_heat_l;
A_= A/B*T_dinf + A/B*nu*Q*Y_0inf/c_p - A*L/c_p;

%surface node
dTdt(noTemperatures) = LiqData.D_heat/delta_r^2*(2*Temperatures(noTemperatures-1) - 2*Temperatures(noTemperatures) - ...
                       (1+1/noTemperatures)*2*delta_r*Temperatures(noTemperatures) +  (1+1/noTemperatures)*2*delta_r*A_);  

%droplet size change
dr2dt = - m_d * 2 * lambda_heat_g / (LiqData.rho *c_p_Air);

%transpose for ode solver
dTdt = dTdt';

Y = [dTdt; dr2dt];

end

