//initialisation Newton SA  + vitesse cible

#ifdef _OPENM4
#pragma omp simd
#endif 

E_Float tol = 1.e-12;
E_Float Cv1 = 7.1;

for (E_Int noind = 0; noind < ifin-ideb; noind++)
{
  // =================================================================================================
  nulocal           = alpha_vec[noind]/yplus_vec[noind]; //recall:yplus_vec[noind ] = rowall*yibc/muwall = yibc/nulocal
  // =================================================================================================
  yibc             = mu_vec[noind]*yplus_vec[noind]/ro_vec[noind];
  yplus            = utau_vec[noind]*yplus_vec[noind];
  yplus_vec[noind] = yplus;
  
  denoml10 = yplus*yplus-8.15*yplus+86.;
  denoml10 = denoml10*denoml10;
  
  grad_linear = 0.;
  if (wl_ibm_swtch==1 || wl_ibm_swtch==11 || wl_ibm_swtch==31){
    umod = utau_vec[noind]*(5.424*atan((2.*yplus-8.15)/16.7) + log10(pow(yplus+10.6,9.6)/denoml10) - 3.52);

    // =================================================================================================
    // Near-wall modification of Spalart-Allmaras turbulence model for IBM, Tamaki, Harada, Imamura, 2017
    // Wall modeling for LES on non-body conforming Cartesian grids, Tamaki & Kawai, 2021
    // =================================================================================================
    // linearization --> Umod(Y)=Umod_image - DUmod/DY|image*(Y_image - Y)  
    // Umod_image = uext_vec[noind]
    // Y_image    = beta_vec[noind]
    // Y          = alpha_vec[noind] 
    // wolfram alpha
    // t = u_tau
    // v = kinematic viscosty
    // y = Y
    // Musker's law u+ = (5.424 ArcTan[(2 y+ - 8.15)/16.7] + Log[(y+ + 10.6)^9.6/((y+)^2 - 8.15 y+ + 86)^2] - 3.52)
    // into wolfram alpha
    // D[t (5.424 ArcTan[(2 t (y/v) - 8.15)/16.7] + Log_10[(t (y/v) + 10.6)^9.6/((t (y/v))^2 - 8.15 t (y/v) + 86)^2] - 3.52), y] 
    // see below for the answer
    // d(u_musker)/dy = t^2 ((7.079 v - 1.73718 t y)/(t^2 y^2 - 8.15 t v y + 86 v^2) + 4.16923/(v ((t y)/v + 10.6)^1) + 0.649581/(v (0.00358564 ((2 t y)/v - 8.15)^2 + 1)))
    // DUmod/DY = utau_vec[noind]**2*(prt1+prt2+prt3)

    prt1        = (7.079*nulocal-1.73718*utau_vec[noind]*beta_vec[noind])/(pow(utau_vec[noind]*beta_vec[noind],2)-8.15*utau_vec[noind]*nulocal*beta_vec[noind]+86*nulocal*nulocal);
    prt2        = 4.16923/(nulocal*(utau_vec[noind]*beta_vec[noind]/nulocal+10.6));
    prt3        = 0.649581/(nulocal*(0.00358564*pow((2*utau_vec[noind]*beta_vec[noind]/nulocal-8.15),2)+1));

    grad_linear = pow(utau_vec[noind],2)*(prt1+prt2+prt3);
  }
  else if (wl_ibm_swtch==2 || wl_ibm_swtch==12 || wl_ibm_swtch==32){
    umod = utau_vec[noind]*(BbarSA_WL + c1SA_WL*log(pow(yplus+a1SA_WL,2)+pow(b1SA_WL,2)) 
			              - c2SA_WL*log(pow(yplus+a2SA_WL,2)+pow(b2SA_WL,2)) 
			              - c3SA_WL*atan2(b1SA_WL,yplus+a1SA_WL)              
			              - c4SA_WL*atan2(b2SA_WL,yplus+a2SA_WL));
    // into wolfram alpha
    // d/dy(t (B + c_1 log(((y t)/v + a_1)^2 + b_1^2) - c_2 log(((y t)/v + a_2)^2 + b_2^2) - c_3 tan^(-1)((y t)/v + a_1, b_1) - c_4 tan^(-1)((y t)/v + a_2, b_2))) = 
    // see below for the answer
    // d(u_sa)/dy = t ((2 c_1 t (a_1 + (t y)/v))/(v ((a_1 + (t y)/v)^2 + b_1^2)) - (2 c_2 t (a_2 + (t y)/v))/(v ((a_2 + (t y)/v)^2 + b_2^2)) + (b_1 c_3 t)/(v ((a_1 + (t y)/v)^2 + b_1^2)) + (b_2 c_4 t)/(v ((a_2 + (t y)/v)^2 + b_2^2)))

    prt1        = (2*c1SA_WL*utau_vec[noind]*(a1SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal))/
                  (nulocal*(pow(a1SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal,2)+pow(b1SA_WL,2)));

    prt2        = (2*c2SA_WL*utau_vec[noind]*(a2SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal))/
                  (nulocal*(pow(a2SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal,2)+pow(b2SA_WL,2)));

    prt3        = (b1SA_WL*c3SA_WL*utau_vec[noind])/
                  (nulocal*(pow(a1SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal,2)+pow(b1SA_WL,2)));

    prt4        = (b2SA_WL*c4SA_WL*utau_vec[noind])/
                  (nulocal*(pow(a2SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal,2)+pow(b2SA_WL,2)));
    
    grad_linear = utau_vec[noind]*(prt1-prt2+prt3+prt4);
  }
  nulocal_vec[noind] = grad_linear;
  umod_lin = uext_vec[noind] - grad_linear*(beta_vec[noind]-alpha_vec[noind]) ;      

  umod     = K_FUNC::E_abs(umod);
  umod     = (1-linearizeWM)*umod+linearizeWM*umod_lin;
  // =================================================================================================

  ucible0 = sign_vec[noind] * umod;
  ucible_vec[noind] += ucible0 * ut_vec[noind]; // vitesse tangentielle pour le pt IBC
  vcible_vec[noind] += ucible0 * vt_vec[noind];
  wcible_vec[noind] += ucible0 * wt_vec[noind];

  // uext: norme de la composante tangentielle de la vitesse externe
  uext = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  uext = std::max(uext, 1.e-12);

  // tcible_vec[noind] = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann
  twall = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);
  tcible_vec[noind] =  twall + (tcible_vec[noind] + 0.5*(uext*uext)/(cv*gamma) - twall)*(umod/uext) - 0.5*(umod*umod)/(cv*gamma); // Equations de Crocco, plus precises pour ecoulements compressibles (Benjamin's formula)

  // van driest pour nut
  expy                = 1.-exp(-yplus/19.);// ranges 17 a 26
  nutcible_vec[noind] = (kappa * alpha_vec[noind])*utau_vec[noind ] * expy*expy;//negatif si pt ibc interieur aux corps
  // printf(" Avant nutcible  pour l indice %d : %g \n", noind, nutcible_vec[noind]);
  E_Float nutcible = K_FUNC::E_abs( nutcible_vec[noind] );
  // equation 4eme degre
  E_Float nu = (mu_vec[noind]/ro_vec[noind]);
  E_Float a = nutcible/nu;
  E_Float b = a*Cv1*Cv1*Cv1;
  //printf("x^4 + %g x^3 +%g = 0\n", -a,-b);

  // debug
  //a = 1.;
  //b = -12.;

  // changement de variable 4eme degre
  E_Float p = -3*a*a/8.;
  E_Float q = -a*a*a/8.;
  E_Float r = -3*a*a*a*a/256.-b;
  //E_Float a_adim = a/(mu_vec[noind]/ro_vec[noind]);

  // equation 3eme degre
  E_Float ap = 8.;
  E_Float bp = -4*p;
  E_Float cp = -8*r;
  //E_Float dp = (a*a)/2.*(pow(a/4.,4.)+3*b);
  E_Float dp = 4*p*r-q*q;
  //printf("cubique: %g x^3 + %g x^2 + %g x + %g\n",ap,bp,cp,dp);

  // racines du 3eme degre
  E_Float delta = 18*ap*bp*cp*dp-4*bp*bp*bp*dp+bp*bp*cp*cp-4*ap*cp*cp*cp-27*ap*ap*dp*dp;
  E_Float delta0 = bp*bp-3*ap*cp;
  E_Float delta1 = 2*bp*bp*bp-9*ap*bp*cp+27*ap*ap*dp;

  E_Float superdelta = -27.*ap*ap*delta;

  //printf("delta %g superdelta>0 = %g\n", delta, superdelta);
  
  
  E_Float y1 = -1.;
  E_Float y2 = -1.;
  if (fabs(delta) < tol && fabs(delta0) < tol) 
  {
    y1 = -bp/(3*ap);  //printf("racine y1=%g\n", y1);
  }
  else if (fabs(delta) < tol)
  {
    y1 = (9*ap*dp - bp*cp)/(2*delta0);
    y2 = (4*ap*bp*cp-9*ap*ap*dp-bp*bp*bp)/(ap*delta0);
    //printf("racine y1=%g y2=%g\n", y1, y2);
  }
  else
  { 
    /* version super delta */
    E_Float C1 = -1.; E_Float C2 = -1.;
    if (superdelta >= 0.)
    {
      E_Float root = sqrt(superdelta);
      if (delta1-root >= 0.)
        C1 = pow( (delta1 -root) /2., 1./3. );
      else
        C1 = -pow( (root -delta1) /2., 1./3. );
      if (delta1 +root >= 0)
        C2 = pow( (delta1 +root) /2., 1./3. );
      else
        C2 = -pow( -(delta1 +root) /2., 1./3. );
    }
  
    //printf("C1=%g, C2=%g\n", C1, C2);
    y1 = -1./(3*ap)*(bp + C1 + delta0/C1 );
    y2 = -1./(3*ap)*(bp + C2 + delta0/C2 );
    //printf("racine y1=%g y2=%g\n", y1, y2);
  }
 

  // racine de l equation du 4eme degre
  E_Float c1 = 2*y1-p;
  E_Float c2 = 2*y2-p;
  //printf("c1 > 0 = %g, c2 > 0 = %g\n",c1,c2);

  E_Float z1 = -123456.;
  if (c1 >= tol)
  {
    E_Float p1 = -2*y1-p+2*q/(sqrt(c1));
    E_Float p2 = -2*y1-p-2*q/(sqrt(c1));
    //printf("1. p1=%g p2=%g\n", p1,p2);
    if (p1 >= tol) { z1 = 0.5*(sqrt(c1)+sqrt(p1));}
    else if (p2 >= tol) { z1 = 0.5*(sqrt(c1)+sqrt(p2));}
  }
  if (c2 >= tol && z1 == -123456)
  {
    E_Float p1 = -2*y2-p+2*q/(sqrt(c2));
    E_Float p2 = -2*y2-p-2*q/(sqrt(c2));
    //printf("2. p1=%g p2=%g\n", p1,p2);
    if (p1 >= tol) { z1 = 0.5*(sqrt(c2)+sqrt(fabs(p1)));}
    else if (p2 >= tol) { z1 = 0.5*(sqrt(c2)+sqrt(fabs(p2)));}
  }
  if (c1 <= tol && z1 == -123456)
  {
    E_Float b0 = y1*y1-r;
    if (b0 >= tol) 
    {
      E_Float p1 = -2*y1-p+4.*sqrt(b0);
      E_Float p2 = -2*y1-p-4.*sqrt(b0);
      //printf("3. p1=%g p2=%g\n", p1,p2);
      if (p1 >= tol) { z1 = 0.5*(sqrt(fabs(c1))+sqrt(fabs(p1))); }
      else if (p2 >= tol) { z1 = 0.5*(sqrt(fabs(c1))+sqrt(fabs(p2))); }
    }
  }
  if (c2 <= tol && z1 == -123456)
  {
    E_Float b0 = y2*y2-r;
    if (b0 >= tol) 
    {
      E_Float p1 = -2*y2-p+4.*sqrt(b0);
      E_Float p2 = -2*y2-p-4.*sqrt(b0);
      //printf("4. p1=%g p2=%g\n", p1,p2);
      if (p1 >= tol) { z1 = 0.5*(sqrt(fabs(c2))+sqrt(fabs(p1))); }
      else if (p2 >= tol) { z1 = 0.5*(sqrt(fabs(c2))+sqrt(fabs(p2))); }
    }
  }

  // nutile final
  E_Float nutilde1 = (z1 + a/4.)*nu;
  aa_vec[noind] = nutilde1;
  // aa_vec[noind] = kappa*utau_vec[noind]*yibc;
  // printf("nutilde final pour l'indice %d = %g\n", noind, nutilde1);

}
