clear;
clc;
format short e;
close all;
fclose all;

decomp_ref = readtable('4.86e+03kgm3_100M_9mu_0.1musec_2300K_referenz.csv');
decomp =  readtable('4.86e+03kgm3_100M_9mu_0.1musec_2300K.csv');

conc_ref = table2array(decomp_ref(5:end, 2));
conc = table2array(decomp(5:end, 2));

figure
hold on
box on
plot(conc_ref);
plot(conc);