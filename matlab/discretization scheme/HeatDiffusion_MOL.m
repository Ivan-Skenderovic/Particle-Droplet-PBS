function [Y] = HeatDiffusion_MOL(t, Temperatures)

D_heat = 1;
T_inf = 400;
h = 1;
DropletRadius = 1;
noTemperatures = length(Temperatures);

%set radial coordinate
r = linspace(0, DropletRadius, noTemperatures);
delta_r = r(2) - r(1);

% symmetry boundary condition at center node
dTdt(1) = D_heat*6/delta_r^2*(-Temperatures(1) + Temperatures(2));

%inner temperature nodes
innerNodeIdx = 2:noTemperatures-1;
dTdt(innerNodeIdx) = D_heat*( ...
(Temperatures(innerNodeIdx-1) - 2.*Temperatures(innerNodeIdx) + Temperatures(innerNodeIdx+1))/delta_r.^2 ...
+ 1./r(innerNodeIdx)'*1./delta_r.*(Temperatures(innerNodeIdx+1) - Temperatures(innerNodeIdx-1)));

%apply boundary condition
dTdt(noTemperatures) = 2*D_heat/delta_r^2*( ...
    Temperatures(noTemperatures-1) - Temperatures(noTemperatures)...
    -(1+1/noTemperatures)*h*delta_r*(Temperatures(noTemperatures) - T_inf) ...
    );

dTdt = dTdt';

Y = dTdt;

end

