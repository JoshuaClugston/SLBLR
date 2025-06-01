# SLBLR
This repository contains a Julia implementation of the Surrogate "Level-Based" Lagrangian Relaxation (SLBLR) method, applied to the generalized assignment problem (GAP). Original SLBLR methodology is implemented (including slight modification) in accordance with

Bragin, M.A., Tucker, E.L. Surrogate “Level-Based” Lagrangian Relaxation for mixed-integer linear programming. 
Sci Rep 12, 22417 (2022). https://doi.org/10.1038/s41598-022-26264-1

using JuMP and CPLEX in Julia. All data files used to test the methododology were acquired using Mutsunori Yakiura’s GAP instances at 
http://www.co.mi.i.nagoya-u.ac.jp/yagiura/gap/, and are included alongside the main scripts of the SLBLR implementation. Particularly, 
the e201600, e801600, and d05100 described therein are used to test the Julia implementation contained in this documentation. Further, a 
JuMP model using CPLEX is provided as a comparison.






