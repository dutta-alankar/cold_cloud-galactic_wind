#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i,j,k;
  double ***T;
  double ***p, ***rho;
  double *dx, *dy;
  double mu = MeanMolecularWeight(d->Vc);
  //printf("%f\n",mu);
  T = GetUserVar("T");
  rho = d->Vc[RHO]; // pointer shortcut to density
  p = d->Vc[PRS];   // pointer shortcut to pressure
  dx = grid->dx[IDIR]; // shortcut to dx 
  dy = grid->dx[JDIR]; // shortcut to dy 
  DOM_LOOP(k,j,i){
    T[k][j][i] = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*KELVIN*mu;
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





