ax = aa_vec[noind]*utau_vec[noind];  
l1 = pow(ax+10.6,9.6);
l2 = ax*(ax-8.15) + 86.;
l3 = (2.*ax-8.15)/16.7;

tp =  aa_vec[noind]*( 9.6*pow((ax + 10.6),8.6)  - l1/l2*(  4.*aa_vec[noind]*utau_vec[noind] - 16.30 ) );

//fp =  5.424*2./16.7*aa/(1. + ((2.*ax - 8.15)/16.7)*((2.*ax - 8.15)/16.7) )+ tp/(l1*log(10.)) + bb/(utau*utau);

fp = 0.649580838323*aa_vec[noind]/(1. + l3*l3 )  +  tp/(l1*2.302585093) + uext_vec[noind]/(utau_vec[noind]*utau_vec[noind]);

// into wolfram alpha
// D[5.424 ArcTan[(2 y t - 8.15)/16.7] + Log[10, (y t + 10.6)^9.6/((y t)^2 - 8.15 y t + 86)^2] - 3.52 - u/t, t]
// answer
// (y (7.079 - 1.73718 t y))/(t^2 y^2 - 8.15 t y + 86) + u/t^2 + (4.16923 y)/(t y + 10.6) + (0.649581 y)/(0.0143426 (t y - 4.075)^2 + 1)
