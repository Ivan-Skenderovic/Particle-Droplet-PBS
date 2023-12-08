function [Y] = HeatDiffusion_1D(t, TandRsquare, LiqData, GasData)

%ConfigSettings

noTemperatures = length(TandRsquare) - 1;
Rsquare = TandRsquare(end);

Temperatures = TandRsquare(1:end-1);

if(sqrt(TandRsquare) > 1e-10)    
    r = sqrt(Rsquare);
    delta_r = r/noTemperatures;
    if ~isreal(delta_r)
        warning('delta_r');
    end
        
    r(1:noTemperatures+1) = 0: delta_r: r;

    dTdt(1:noTemperatures) = 0;
    innerNodeIdx = 2:noTemperatures - 1;

    % center temperature node
    dTdt(1) = LiqData.D_heat*6/delta_r^2*(-Temperatures(1) + Temperatures(2));

    %inner temperature nodes
    dTdt(innerNodeIdx) = LiqData.D_heat*( ...
    (Temperatures(innerNodeIdx-1) - 2.*Temperatures(innerNodeIdx) + Temperatures(innerNodeIdx+1))/delta_r.^2 ...
    + 1./r(innerNodeIdx)'*1./delta_r.*(Temperatures(innerNodeIdx+1) - Temperatures(innerNodeIdx-1)));

    surfaceTemperature = Temperatures(end);
        
    calcTransportProperties
    
    B = (T_dinf - T_ds + nu*Y_0inf*Q/L)/H_d;        

    A = log(1+B)./r(noTemperatures)*lambda_heat_g/lambda_heat_l;
    
    A_= A/B*T_dinf + A/B*nu*Q*Y_0inf/c_p - A*L/c_p;

    dTdt(noTemperatures) = LiqData.D_heat/delta_r^2*(2*Temperatures(noTemperatures-1) - 2*Temperatures(noTemperatures) - ...
                           (1+1/noTemperatures)*2*delta_r*Temperatures(noTemperatures) +  (1+1/noTemperatures)*2*delta_r*A_);  
   
    dr2dt = - m_d * 2 * lambda_heat_g / (LiqData.rho *c_p_Air);

    dTdt = dTdt';

    Y = [dTdt; dr2dt];
else

dTdt(1:noTemperatures) = 0;
dr2dt = 0;

Y = [dTdt'; dr2dt];    

end

