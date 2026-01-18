(* This is a mathematica script to generate the viscous
tensor coefficients for the Bjorken 2D expansion case 

The output is the C++ files:
  knx_shear.txt
  kxx_shear.txt
  knn_shear.txt
  knx_bulk.txt
  kxx_bulk.txt
  knn_bulk.txt

These files are procesed by the processbj2d.cpp code to generate VischydroNode_inc.hpp

Variables:
  vx, vy : fluid velocity components in x and y directions
  v2 = vx^2 + vy^2
  cs2 : speed of sound squared

Author: D. Teaney
Date: Jan 10, 2024
*)
dij = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1} }

vv = { vx, vy, 0 } 
hij = Simplify[ dij - Outer[ Times, vv, vv ]  ]

g = 1/Sqrt[1 - (vx^2 + vy^2) ]

hh =  Outer[ Times, hij, hij  ]  

himjn  = Simplify[(Transpose[hh, 2<->3 ] + Transpose[hh, 1<->3])/2 ]

uiuj  = Simplify[ Outer[ Times, vv, vv ]  g^2  ]
dij = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1} }

P = FullSimplify[ -cs2 uiuj + 1/3 (dij + uiuj)   ] 

hPh = Simplify[ hij.P.hij ]  

x1 = Simplify[ Tr[P.hij]  ]  
x2 = Simplify[ Tr[P.hij.P.hij ]  ]  
x3 = Simplify[x2/x1^2]

t1 = Simplify[ Outer[Times, hij, hPh]/x1 ]  
t2 = Simplify[ Outer[Times, hPh, hij]/x1 ]  

d=3
shear =  FullSimplify[ FullSimplify[ 2 (himjn - t1  - t2 + x3 Simplify[hh]) ]   /. {vx^2 + vy^2 -> v2 }  ]

knx = Export["knx_shear.txt", Flatten[ FullSimplify[(1- cs2 v2)^2 shear[[1;;2, 1;;2, 3, 3]] ] ] // CForm] 
kxx = Export["kxx_shear.txt", Flatten[ FullSimplify[(1- cs2 v2)^2 shear[[1;;2, 1;;2, 1;;2, 1;;2]] ] ] // CForm] 
knn = Export["knn_shear.txt", FullSimplify[(1- cs2 v2)^2 shear[[3, 3, 3, 3]]  ] // CForm]

bulk = FullSimplify[ FullSimplify[ hh/x1^2] /. {vx^2 + vy^2 -> v2 } ]

knxB = Export["knx_bulk.txt", Flatten[ FullSimplify[(1- cs2 v2)^2 bulk[[1;;2, 1;;2, 3, 3]] ] ] // CForm] 
kxxB = Export["kxx_bulk.txt", Flatten[ FullSimplify[(1- cs2 v2)^2 bulk[[1;;2, 1;;2, 1;;2, 1;;2]] ] ] // CForm] 
knnB = Export["knn_bulk.txt", FullSimplify[(1- cs2 v2)^2 bulk[[3, 3, 3, 3]]  ] // CForm]

