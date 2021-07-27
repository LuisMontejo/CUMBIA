# CUMBIA
CUMBIA - a set matlab/octave codes for the analysis of reinforced concrete members

cite as: Montejo, L. A., & Kowalsky, M. J. (2007). CUMBIA—Set of codes for the analysis of reinforced concrete members. CFL technical rep. no. IS-07, 1.

## Change log for v0.2

Calculation of Mn was revised (there was a problem for some sections where Mn was not calculated in v0.1)

Calculation of Mi(k) in the axial-moment interaction was revised (there was a problem for some sections where Mi(k) was not calculated in v0.1)

Axial force values for the calculation of the axial-moment interaction curve were changed (depending on the section details, for compression axial load ratios greater than 70% the algorithm was not converging)

The estimation of first yield curvature now includes the compression strain in the concrete (1.8*fpc/Ec), before it only checked for tensile strain in the rebar (fy/Es)
“breaks” outside a loop were replaced with “returns” 
