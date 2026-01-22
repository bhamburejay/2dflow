(*   Code to convert from Density to Landau Frame 

The purpose here is to generate some code for converting between the Landau
and density Frames.  We are solving the landau frame condition

The code assumes that the stress Ttt, Ttx, Tty, Txx, Txy, Tyy are defined
by the user, and that ux and uy are the fluid velocities in the x and y directions respectively.

It then computes the function Fx and Fy 
and their derivatives Mxx, Mxy, Myx, Myy

Mxy = D[Fx, uy] and Myx = D[Fy, ux]

The outputs are written to the files DFtoLF_Func_inc.hpp and DFtoLF_DFunc_inc.hpp
for the function and its derivatives respectively.

These functions contain ccode 

Fx=...;
Fy=...;

Mxx=...;
Mxy=...;
Myx=...;
Myy=...;

Author: Derek Teaney
Date: January 2026

*)

ut = - Sqrt[1 + ux ^ 2 + uy ^ 2]; 


eg = Ttt ut ^ 2 + 2 Ttx ut ux + 2 Tty ut uy + Txx ux ^ 2 + 2 Txy ux uy + Tyy uy ^ 2;

Fx = Ttx ut + Txx ux + Txy uy + eg ux
Fy = Tty ut + Txy ux + Tyy uy + eg uy 

Print[ Fx, Fy]
m = FullSimplify[ D[ {Fx, Fy}, {{ux, uy}} ] ] 
(* Simplify[ D[ {Fx, Fy}, ux, uy ]  ] *)

FullSimplify[m[[2, 1]]- D[ Fy, ux] ]



<<Format.m 

(* Open a file  DFtoLF_inc.hpp *)
SetOptions[CAssign, AssignOptimize->False, AssignBreak->False]
stream = OpenWrite["DFtoLF_Func.inc"]
(* 
Write the assignments to the file
*)
(* SetOptions[Print, PageWidth -> Infinity, OutputStream -> stream]  *)

Write[stream, CAssign["Fx", Fx ]  ]
Write[stream, CAssign["Fy", Fy ]  ]
Close[stream]

stream = OpenWrite["DFtoLF_Jac.inc"]

Write[stream, CAssign["Mxx", m[[1, 1]] ] ] 
Write[stream, CAssign["Mxy", m[[1, 2]] ] ]
Write[stream, CAssign["Myx", m[[2, 1]] ] ]
Write[stream, CAssign["Myy", m[[2, 2]] ] ]
Close[stream] 
