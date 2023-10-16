function [value,isterminal,direction] = dropletEvents(t,y)
    %watch droplet diameter
    value = y(end) - 1e-16; 
    isterminal = 1;
    direction = 0;
end

