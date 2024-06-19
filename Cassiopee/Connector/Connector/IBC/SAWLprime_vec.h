ax = aa_vec[noind]*utau_vec[noind];  

prt1SAWL_prime= 2.*c1SA_WL*aa_vec[noind]*(a1SA_WL+ax)/(pow(a1SA_WL+ax,2)+pow(b1SA_WL,2));
prt2SAWL_prime= 2.*c2SA_WL*aa_vec[noind]*(a2SA_WL+ax)/(pow(a2SA_WL+ax,2)+pow(b2SA_WL,2));
prt3SAWL_prime= b1SA_WL*c3SA_WL*aa_vec[noind]/(pow(a1SA_WL+ax,2)+pow(b1SA_WL,2));
prt4SAWL_prime= b2SA_WL*c4SA_WL*aa_vec[noind]/(pow(a2SA_WL+ax,2)+pow(b2SA_WL,2));
fp           = prt1SAWL_prime -
               prt2SAWL_prime +
               prt3SAWL_prime +
               prt4SAWL_prime +
               uext_vec[noind]/(utau_vec[noind]*utau_vec[noind]);

// into wolfram alpha: Log = ln
// d/dt(B + c_1 log((y t + a_1)^2 + b_1^2) - c_2 log((y t + a_2)^2 + b_2^2) - c_3 tan^(-1)(y t + a_1, b_1) - c_4 tan^(-1)(y t + a_2, b_2) - u/t)
// answer::
// (2 c_1 y (a_1 + t y))/((a_1 + t y)^2 + b_1^2) - (2 c_2 y (a_2 + t y))/((a_2 + t y)^2 + b_2^2) + (b_1 c_3 y)/((a_1 + t y)^2 + b_1^2) + (b_2 c_4 y)/((a_2 + t y)^2 + b_2^2) + u/t^2

// from Allmaras et al. 2012 paper:
// u+(y+) = B + c1 log(y+ + a1)2 + b2 1 − c2 log (y+ + a2)2 + b2 ? ? 2 − c3ArcTan[y+ + a1, b1] − c4ArcTan[y+ + a2, b2]
// where ArcTan[x, y] is the Mathematica function equivalent atan2 fortran
// Mathematic and Wolfram Alpha have the same parent company and both use the wolfram language.


