ax = aa_vec[noind]*utau_vec[noind]; 

utauv_vec[noind] = BbarSA_WL + c1SA_WL*log(pow(ax+a1SA_WL,2)+pow(b1SA_WL,2)) 
                             - c2SA_WL*log(pow(ax+a2SA_WL,2)+pow(b2SA_WL,2)) 
                             - c3SA_WL*atan2(b1SA_WL,ax+a1SA_WL)              
                             - c4SA_WL*atan2(b2SA_WL,ax+a2SA_WL) - uext_vec[noind]/utau_vec[noind];

if (K_FUNC::E_abs( utauv_vec[noind] ) > newtoneps && skip == 0) err = 0;
//see Modificaations and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model
//to better understand the atan2 arguments order... Allmaras writes ArcTan[x,y] but c++ needs atan2(y,x)
//Allmaras clarifies that ArcTan[x,y] is the Mathematica equivalent of atan2(y,x)
//therefore when Allmaras writes ArcTan[yplus+a1,b1], x = yplus+a1 and y = b1 and the fortran function is
//atan2(b1,yplus+a1)
//Note: checked through a python script the correct implemenation


