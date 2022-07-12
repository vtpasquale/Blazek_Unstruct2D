# Blazek Unstruct2D

This is a modified version of `Unstruct2D` by Blazek [1]. Code changes from the original version are tracked with Git. Feature change(s):
* Volume output data format changed from Vis2D format to legacy VTK format. The new format can be read using Paraview.
* Select volume output [rho,u,v,w,p] is also written to a Gamma-format [2] solution file (.sol or solb). 
  - This output is intended for use with refine [3]. 
  - For refine compatibility, Gamma-format output uses FUN3D nondimensionalization.

[1] Blazek, Jiri. _Computational fluid dynamics: principles and applications_. Butterworth-Heinemann, 2015.
[2] https://github.com/LoicMarechal/libMeshb
[3] https://github.com/nasa/refine