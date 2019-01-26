/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "math.h"
#include "coolit.h"
#include <stdbool.h>
//#include "gsl.h"

double mu, Rcl, tcc, mcl, tdrag, vwind;
int counter = 0;
//velunit kms^-1    dist unit pc    time unit 3.086e+13s
/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  double x0 = g_inputParam[X0];
  double y0 = g_inputParam[Y0];
  double z0 = g_inputParam[Z0];
  double cs;
  double chi = g_inputParam[CHI];
  double Tcl = g_inputParam[TCL];
  double mach = g_inputParam[MACH];
  double alpha = g_inputParam[ALPHA];
  double Ncl = g_inputParam[NCL];
  
  double Twind = chi*Tcl;
  double Tmix = sqrt(Tcl*Twind);
  mu = MeanMolecularWeight(v);
  double P3 = (Ncl/1.0e3)*Tcl;
  double lambdamix = lambda(Tmix)/pow(10,-21.4);
  
  Rcl = 2*(pow(Tcl/1.0e4,(5./2.))*mach/(P3*lambdamix))*(chi/100)*pow(alpha,-1); //pc
  mcl = ((4*CONST_PI/3)*pow(Rcl*UNIT_LENGTH,3.0)*(Ncl*mu*UNIT_DENSITY))/CONST_Msun; //solarmass
  
  v[TRC] = 0.0;
  v[PRS] = Ncl*mu*Tcl/(KELVIN*mu);
  v[VX2] = 0.;
  v[VX3] = 0.; 
  double r;
  r = sqrt((x1-x0)*(x1-x0) + (x2-y0)*(x2-y0) + (x3-z0)*(x3-z0));
  if (r <= Rcl){ //cold cloud
    v[RHO] = Ncl*mu;
    v[VX1] = 0.;
    //v[PRS] = v[RHO]*Tcl/(KELVIN *mu);
    v[TRC] = 1.0;
    //cs = sqrt(g_gamma*v[PRS]/v[RHO]);
  }
  else{
    v[RHO] = Ncl*mu/chi;
    //v[PRS] = v[RHO]*Twind/(KELVIN *mu);
    cs = sqrt(g_gamma*v[PRS]/v[RHO]);
    v[VX1] = mach*cs;
    vwind = mach*cs;
  }
  //#if HAVE_ENERGY
  //#endif
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  tcc = (pow(chi,0.5)*Rcl*(CONST_pc*1e-5)/(cs*mach))/(UNIT_LENGTH/UNIT_VELOCITY); //code units
  tdrag = sqrt(chi)*tcc; //code units
  if (counter==0 && world_rank==0) {
    printf("Cloud size is %e pc and Cloud Crushing time is %e s (%e [code units]) \nInitial Pressure set to %e [code units]\n",Rcl,pow(chi,0.5)*Rcl*(CONST_pc*1e-5)/(cs*mach),(pow(chi,0.5)*Rcl*(CONST_pc*1e-5)/(cs*mach))/(UNIT_LENGTH/UNIT_VELOCITY),Ncl*mu*Tcl/(KELVIN*mu));
    //printf("Other useful stuff: \nKELVIN: %e, \n",);
    counter++;
  }
  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

  int i, j, k;
  int count_hot = 0;
  int count_cold = 0;
  double Tcl = g_inputParam[TCL];
  double Ncl = g_inputParam[NCL];
  double chi = g_inputParam[CHI];
  double Twind = chi*Tcl;
  double Tmix = sqrt(Tcl*Twind);
  double dV, vol, Temperature;
  double mass, vx2, vy2, vz2, speed, speed_hot = 0., speed_cold = 0., dragf, mass_all, speed_mix = 0., speed_mix_all;
  double *dx, *dy, *dz;

  /* ---- Set pointer shortcuts ---- */
  dx = grid->dx[IDIR];
  dy = grid->dx[JDIR];
  dz = grid->dx[KDIR];
  
  /* ---- Main loop ---- */
  mass = 0.0; speed = 0.0;
  DOM_LOOP(k,j,i){
    dV = dx[i]*dy[j]*dz[k]; // Cell volume (Cartesian coordinates)
    if (d->Vc[RHO][k][j][i]>Ncl*mu/3.0) mass += d->Vc[RHO][k][j][i]*dV; 

    vx2 = d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i];
    vy2 = d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i];
    vz2 = d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i];
    speed = sqrt(vx2+vy2+vz2);
    Temperature = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*KELVIN*mu;
    //printf("%e\n",Temperature);
    if (Temperature<=Tmix)
    {
      if (speed>speed_mix) speed_mix = speed;
    }
    
  }
  mass = mass*UNIT_DENSITY*pow(UNIT_LENGTH,3.0)/CONST_Msun;
  //if(count_cold!=0) speed_cold = speed_cold/count_cold;
  //if(count_hot!=0) speed_hot = speed_hot/count_hot;

  /* ---- Parallel data reduction ---- */
  #ifdef PARALLEL
  int size; 
  //int world_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  //MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Allreduce (&mass, &mass_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&speed_mix, &speed_mix_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  //MPI_Allreduce (&speed_hot, &speed_hot_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  #else
  speed_mix_all = speed_mix;
  mass_all = mass;
  #endif
  mass_all /= mcl;
  dragf = (vwind-speed_mix_all)/vwind;
  if (mass_all<0.1) return;

  /* ---- Write ascii file "analysis.dat" to disk ---- */
  if (prank == 0)
  {
    char fname[512];
    static double tpos = -1.0;
    FILE *fp;
    sprintf (fname, "%s/analysis.dat",RuntimeGet()->output_dir);
    if (g_stepNumber == 0)
    { // Open for writing only when we’re starting
      fp = fopen(fname,"w"); /* from beginning */
      fprintf (fp,"Chi: %f\ttcc [code units]: %e\tMcl [Solar mass]: %e\tmu: %f\t [code time] = %e s\n\n", chi, tcc, mcl,mu, UNIT_LENGTH/UNIT_VELOCITY);
      fprintf (fp,"%s\t\t%s\t\t%s\t%s\n", "t/tcc", "t/tdrag", "M(rho>rhocl/3)/Mcl","deltav/vwind");
    }
    else
    {
      /* Append if this is not step 0 */
      if (tpos < 0.0)
      {      // Obtain time coordinate of to last written row
        char sline[512];
        fp = fopen(fname,"r");
        while (fgets(sline, 512, fp)) {}
        sscanf(sline, "%lf\n",&tpos); // tpos = time of the last written row 
        fclose(fp);
      }
      fp = fopen(fname,"a");
  }
  if (g_time > tpos)
    { // Write if current time if > tpos
      fprintf (fp, "%e\t%e\t%e\t\t%e\n",g_time/tcc, g_time/tdrag,mass_all, dragf);
    }
    fclose(fp);
  }
}


#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double cs;
  double Tcl = g_inputParam[TCL];
  if (side == 0){
    TOT_LOOP(k,j,i){
      if ((d->Vc[PRS][k][j][i]*KELVIN *mu/d->Vc[RHO][k][j][i]) < Tcl) {  //artificial heating (Temperature floor)
        d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*Tcl/(KELVIN *mu);
      }
    }
  }

  //x1 = grid->x[IDIR];  /* -- array pointer to x1 coordinate -- */
  //x2 = grid->x[JDIR];  /* -- array pointer to x2 coordinate -- */
  //x3 = grid->x[KDIR];  /* -- array pointer to x3 coordinate -- */

  if (side == X1_BEG){    /* -- select the boundary side -- */
    BOX_LOOP(box,k,j,i){  /* -- Loop over boundary zones -- */
      d->Vc[RHO][k][j][i] = g_inputParam[NCL]*mu/g_inputParam[CHI];
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*g_inputParam[CHI]*g_inputParam[TCL]/(KELVIN *mu);
      cs = sqrt(g_gamma*d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]);
      d->Vc[VX1][k][j][i] = g_inputParam[MACH]*cs;
      d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[VX3][k][j][i] = 0.0;
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
