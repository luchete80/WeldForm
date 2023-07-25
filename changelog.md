# Changelog

## [0.4.1] - TBD

### Solver 
  - Now reading constant thermal properties and initial conditions. First thermal-mechanical example ready.
  - Fixed Friction Stir Welding example (greater particle density results in not so negative temp at the begining)
  - Added input reading details (to know in which vars exists an input error)
  - Randles and Libersky variable time step integrator working
### Input File 
  - Nb frequency update option on input 
  - 
  - Added Cylindrical grid option in input 
  ### Deprecated
 - Original Simultaneous Solver (All vars calculated at same time)

## [0.4.0] - 2023-02-03
### Core 
 - Structural , thermal and coupled solver.
 - Fraser, Leapfrog, Verlet and Differential KickDrift time integration working.
 - Bilinear, Hollomon and Johnson Cook materials.
 - Frictional Wang contact.
 - Working accumulated frictional contact Work, external forces Work and plastic dissipated Work.  
 
### Input File 
 - Wang Frictional contact algorithm working, giving same results that FEM solvers.
 - First imput file reader: working only for structural problems.
 - Input Boundary conditions of velocity and zones. 
 


