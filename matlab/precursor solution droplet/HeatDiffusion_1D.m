function dYdt = HeatDiffusion_1D(t, TandRsquare, LiqData, GasData)

noTemperatures = length(TandRsquare) - 1;
Rsquare = TandRsquare(end);
Temperatures = TandRsquare(1:end-1);
DropletRadius = sqrt(Rsquare);
if DropletRadius < 1e-9 %m, termination event at lower droplet size limit
    dYdt(1:length(TandRsquare), 1) = 0;
    return
end
%discretize radial coordinate
r = linspace(0, DropletRadius, noTemperatures);
delta_r = r(2) - r(1);

%center node 
dTdt(1) = LiqData.D_heat*6/delta_r^2*(-Temperatures(1) + Temperatures(2));
innerNodeIdx = 2:noTemperatures - 1;
%inner nodes
dTdt(innerNodeIdx) = LiqData.D_heat*( ...
(Temperatures(innerNodeIdx-1) - 2.*Temperatures(innerNodeIdx) + Temperatures(innerNodeIdx+1))/delta_r.^2 ...
+ 1./r(innerNodeIdx)'*1./delta_r.*(Temperatures(innerNodeIdx+1) - Temperatures(innerNodeIdx-1)));
       
%specific gas constant
Rg = LiqData.Rg/LiqData.calgToJkg; %cal/gm
%specific latent heat of vaporization
L = LiqData.L/LiqData.calgToJkg; %cal/gm
%heat of reaction per unit mass of fuel reacted
Q = LiqData.Q/LiqData.calgToJkg; %cal/gm
%stochiometric fuel oxidizer mass ratio
nu = LiqData.nu;
T_boil = LiqData.T_boil; %K
p_inf = GasData.p_inf; %
%oxidizer mass fraction
Y_0inf = GasData.Y_0inf;
YFf = 0; % assuming complete combustion
T_dinf = GasData.Cg_Law/L*GasData.T_inf; 
surfaceTemperature = Temperatures(end);
T_ds = GasData.Cg_Law/L.*surfaceTemperature;
T_db = GasData.Cg_Law/L.*T_boil;
M_air = GasData.M_air; %kg/mol
W_A = M_air; %average of all species except fuel in kg/mole
W_F = LiqData.M; %kg/mol

%Fuel vapor fraction at droplet surface
Y_Fs = (1 + (W_A/W_F)*(p_inf*exp((GasData.Cg_Law/Rg)*(1./T_ds - 1./T_db))-1)).^-1;      
H_d = (1 - Y_Fs).*(T_dinf - T_ds + nu.*Y_0inf.*Q./L) ./ (Y_Fs + nu*Y_0inf - YFf.*(1 + nu*Y_0inf));
%dimensionless evaporation rate   
m_d = log(1 + (T_dinf - T_ds + nu.*Y_0inf.*Q./L)./H_d);     
% Spalding heat transfer number
B = (T_dinf - T_ds + nu*Y_0inf*Q/L)/H_d;
%Abbreviation variable 
A = log(1+B)./r(noTemperatures)*GasData.lambda_heat/LiqData.lambda_heat;
%Abbreviation variable
A_= A/B*T_dinf + A/B*nu*Q*Y_0inf/GasData.Cg_Law - A*L/GasData.Cg_Law;

%surface node
dTdt(noTemperatures) = LiqData.D_heat/delta_r^2*(2*Temperatures(noTemperatures-1) - 2*Temperatures(noTemperatures) - ...
                       (1+1/noTemperatures)*2*delta_r*Temperatures(noTemperatures) + (1+1/noTemperatures)*2*delta_r*A_);  

dRadius2dt = - m_d * 2 * GasData.lambda_heat / (LiqData.rho*GasData.C_g);

dYdt = [dTdt'; dRadius2dt];    

end

