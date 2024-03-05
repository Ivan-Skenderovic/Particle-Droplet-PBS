% convert volume fractions to mole fractions in an ideal liquid mixture
clear;
clc;

rho_EtOH = 789.45; % kg/m³, density
M_EtOH = 0.046069; % kg/mole, molar mass

rho_2EHA = 903; % kg/m³, density
M_2EHA = 0.144214; % kg/mole, molar mass

% density = mass / volume, volume = mass / density
% mass = M * mole
% volume = M / density 
n_EtOH = 1;
n_2EHA = 0.679; 
x_EtOH = n_EtOH / (n_EtOH + n_2EHA);
%x_2EHA = 1 - x_EtOH;

v_EtOH =  ( n_EtOH * M_EtOH/rho_EtOH ) /...
    ( n_EtOH * M_EtOH/rho_EtOH + n_2EHA * M_2EHA/rho_2EHA );