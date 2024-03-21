Particle-Droplet-PBS (PD-PBS) is a Matlab based software to simulate particle accumulation -and growth inside spray droplets. Particle growth is calculate by solving 
the population-balance equations for coagulation and nucleation. Two standalone case studies are provided to demonstrate possible applications in the field of 
particle synthesis and fuel research. It implements the particle models to estimate time of droplet breakup described in "A population-balance 
method for simulation of particle-induced droplet breakup in spray flame synthesis and suspension spray combustion" by Skenderovic et al., 2023. 
 
The software can be run be either executing the 'main.m' script from the Matlab GUI in the precursor_solution_droplet or the suspension_droplet folder or by 
executing the 'run.sh' script using bash. First the user needs to edit the 'ConfigSettings.m' files to implement different materials or to set the
required accuracy. The simulation results are provided in the form of a particle size distribution, among other statistics, in a 'xxx.csv' file. 
Its' filename contains material parameters set by the use. Finally, the output file is stored inside the 'results' folder. 

The user is encouraged to further tailor the code to their needs. In droplet combustion research models with varying degrees of accuracy and computation time
are used depending on the requirements of the users. In this case, the custom models provide droplet size change and temperature data to PD-PBS. Both case 
studies demonstrate how the interface with custom models can be implemented. 

 

