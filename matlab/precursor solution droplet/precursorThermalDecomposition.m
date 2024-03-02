function concentrationChange = precursorThermalDecomposition(concentration, ...
    decompositionRate, timestep)

% Analytical solution for the equation: dC/dt = -decompositionRate. A
% straighforward extension is to define decompositionRate as function of
% temperature.
newConcentration = concentration.*exp(-decompositionRate.*timestep);
concentrationChange = concentration - newConcentration;

end