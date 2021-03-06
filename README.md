# 3Dtumorsimulator
We simulate the spatiotemporal dynamics of a multicellular tumor using a three-dimensional 
(3D) agent-based model, where cells are represented as interacting objects. 
Cells grow stochastically as a function of time and divide upon reaching a critical size. 
The cell-to-cell interaction is characterized by repulsive and adhesive forces. 
Cell-to-cell and cell-to-matrix damping to account for 
the effects of friction experienced by a moving cell because of other cells and by the 
extracellular matrix (ECM) (or collagen matrix), respectively, are considered.  
Each cell is modeled as a deformable sphere with a time-dependent radius. 
Several physical properties, such as the radius, elastic modulus, 
membrane receptor and ligand concentration, and adhesive interaction, characterize each cell.
For more details see: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.021025


%%%%%%%%%

main.m is main code implementing cell growth, dynamics, and birth-death processes \

Force.m is the code which at each time step calculates the force experienced by cells \ 

initCubicGrid.m specifies the initial position of cells \

randgaussrad.m is generates normally distributed random numbers shifted by the mean and standard deviation

%%%%%%%%%%%

Installation

This code was written to run on Matlab version R2018a

Download all the .m files and place them in the same folder. 
