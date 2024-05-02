//GMSH project
// 
//           7
//          | 
//        5 | 6
//  ________/
// / 3    4 
// 1 2
// |
// 0

lc = 1e-3;
thck = 3.0e-4;
ysal = 0.01;
rsal = 0.002;
rint = 0.002;

Point(1)  = {0.0,-ysal,0.0,lc}; 
Point(2)  = {0.0,-rsal,0.0,lc}; 
Point(3)  = {rsal,-rsal,0.0,lc}; 
