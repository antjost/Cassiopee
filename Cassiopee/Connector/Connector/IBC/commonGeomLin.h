// PW: point paroi
// PI: point interpole (ou exterieur)
// IN: (u,v,w) au point interpole
b0 = xPI[noind+ideb]-xPW[noind+ideb];
b1 = yPI[noind+ideb]-yPW[noind+ideb];
b2 = zPI[noind+ideb]-zPW[noind+ideb];

normb = sqrt(b0*b0+b1*b1+b2*b2);
normb = std::max(normb, 1.e-12);

//@ Image Pnt
n0 = b0/normb;
n1 = b1/normb;
n2 = b2/normb;

// uscaln: u ⋅ n
u = uOut[indR];
v = vOut[indR];
w = wOut[indR];

uscaln = u*n0 + v*n1 + w*n2;
  
//composante normale de la vitesse
un = uscaln*n0;
vn = uscaln*n1;
wn = uscaln*n2;

//composante tangentielle de la vitesse 
ut = u-un;
vt = v-vn;
wt = w-wn;

// uext: norme de la composante tangentielle de la vitesse 
uext = sqrt(ut*ut+vt*vt+wt*wt);

t0=ut/uext;
t1=vt/uext;
t2=wt/uext;
