if (wl_ibm_swtch > 30)
  {
    fsLin       = max((beta_vec[noind]-alpha_vec[noind])/beta_vec[noind],0.);
    coefLin     = roOut[indR]*pow(0.41*beta_vec[noind],2);
    grad_linear = nulocal_vec[noind];    
#include "commonGeomLin.h"
    //t_i,j                                       (ti nj + tj ni)
    t11Ptr[noind+ideb]=coefLin*grad_linear*abs(grad_linear)*fsLin*(t0*n0+t0*n0);
    t12Ptr[noind+ideb]=coefLin*grad_linear*abs(grad_linear)*fsLin*(t0*n1+t1*n0);
    t22Ptr[noind+ideb]=coefLin*grad_linear*abs(grad_linear)*fsLin*(t1*n1+t1*n1);
    t13Ptr[noind+ideb]=coefLin*grad_linear*abs(grad_linear)*fsLin*(t0*n2+t2*n0);
    t23Ptr[noind+ideb]=coefLin*grad_linear*abs(grad_linear)*fsLin*(t1*n2+t2*n1);
    t33Ptr[noind+ideb]=coefLin*grad_linear*abs(grad_linear)*fsLin*(t2*n2+t2*n2);
    if (isSkipFlowFieldWMLES==0){
      t11m[indR]=t11Ptr[noind+ideb];
      t12m[indR]=t12Ptr[noind+ideb];
      t22m[indR]=t22Ptr[noind+ideb];
      t13m[indR]=t13Ptr[noind+ideb];
      t23m[indR]=t23Ptr[noind+ideb];
      t33m[indR]=t33Ptr[noind+ideb];     
    }   
  }



//if (wl_ibm_swtch == 31)
//{
//prt1        = (7.079*nulocal_vec[noind]-1.73718*utau_vec[noind]*beta_vec[noind])/(pow(utau_vec[noind],2)*pow(beta_vec[noind],2)-8.15*utau_vec[noind]*nulocal_vec[noind]*beta_vec[noind]+86*nulocal_vec[noind]*nulocal_vec[noind]);
//prt2        = 4.16923/(nulocal_vec[noind]*(utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind]+10.6));
//prt3        = 0.649581/(nulocal_vec[noind]*(0.00358564*pow((2*utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind]-8.15),2)+1));
//grad_linear = pow(utau_vec[noind],2)*(prt1+prt2+prt3);
//}
//else if (wl_ibm_swtch == 32)
//{
//prt1        = (2*c1SA_WL*utau_vec[noind]*(a1SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind]))/
//(nulocal_vec[noind]*(pow(a1SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind],2)+pow(b1SA_WL,2)));
//
//prt2        = (2*c2SA_WL*utau_vec[noind]*(a2SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind]))/
//(nulocal_vec[noind]*(pow(a2SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind],2)+pow(b2SA_WL,2)));
//
//prt3        = (b1SA_WL*c3SA_WL*utau_vec[noind])/
//(nulocal_vec[noind]*(pow(a1SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind],2)+pow(b1SA_WL,2)));
//
//prt4        = (b2SA_WL*c4SA_WL*utau_vec[noind])/
//(nulocal_vec[noind]*(pow(a2SA_WL+utau_vec[noind]*beta_vec[noind]/nulocal_vec[noind],2)+pow(b2SA_WL,2)));
//
//grad_linear = utau_vec[noind]*(prt1-prt2+prt3+prt4);
//}

