if (linearizeWMLES==1)
  {
    fsLin       = max((beta_vec[noind]-alpha_vec[noind])/beta_vec[noind],0.);
    coefLin     = roOut[indR]*pow(0.41*beta_vec[noind],2);
    grad_linear = nulocal_vec[noind];    
#include "commonGeomLin.h"
    //t_i,j                                       (ti nj + tj ni)
    t11Ptr[noind+ideb]=110000+noind+ideb;//coefLin*grad_linear*abs(grad_linear)*fsLin*(t0*n0+t0*n0);
    t12Ptr[noind+ideb]=120000+noind+ideb;//coefLin*grad_linear*abs(grad_linear)*fsLin*(t0*n1+t1*n0);
    t22Ptr[noind+ideb]=220000+noind+ideb;//coefLin*grad_linear*abs(grad_linear)*fsLin*(t1*n1+t1*n1);
    t13Ptr[noind+ideb]=130000+noind+ideb;//coefLin*grad_linear*abs(grad_linear)*fsLin*(t0*n2+t2*n0);
    t23Ptr[noind+ideb]=230000+noind+ideb;//coefLin*grad_linear*abs(grad_linear)*fsLin*(t1*n2+t2*n1);
    t33Ptr[noind+ideb]=330000+noind+ideb;//coefLin*grad_linear*abs(grad_linear)*fsLin*(t2*n2+t2*n2);
  }
