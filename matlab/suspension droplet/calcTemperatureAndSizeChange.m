odesystemhandle = @(t,Y) HeatDiffusion_1D(t, Y, LiqData, GasData);

[times, TemperaturesAndR2] = ode15s(odesystemhandle, [0 TIMESTEP], ...
    initialTemperaturesAndR2, options);

initialTemperaturesAndR2 = TemperaturesAndR2(end,:); %updated for next iteration
%calculate volume weighted mean temperature

dr2dt = TemperaturesAndR2(end); %get regression rate [m²/s]
r_new = sqrt(r^2 - dr2dt*TIMESTEP) 
