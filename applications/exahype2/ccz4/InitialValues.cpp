#include "InitialValues.h"
#include "Properties.h"
/**
 *  *  * This file is automatically created by Peano. I need it to interact with
 *   *   * the Python API, i.e. to read out data set there.
 *    *    */
#include "Constants.h"

#include "tarch/la/Vector.h"

#include "peano4/utils/Globals.h"

#include <stdio.h>
#include <string.h>

#ifdef IncludeTwoPunctures
#include "libtwopunctures/TP_bindding.h"
#endif

void applications::exahype2::ccz4::gaugeWave(
  double * __restrict__ Q, // Q[64+0],
  const tarch::la::Vector<Dimensions,double>&  x,
  double t
) {
  constexpr int nVars = 59;
  constexpr double pi = M_PI;
  constexpr double peak_number = 2.0;
  constexpr double ICA = 0.1; ///< Amplitude of the wave 
  double HH     = 1.0 - ICA*sin( peak_number*pi*( x[0] - t));
  double dxH    = -peak_number*pi*ICA*cos( peak_number * pi*(x[0] - t));
  double dxphi  = - pow(HH,(-7.0/6.0))*dxH/6.0;
  double phi    = pow(( 1.0 / HH),(1.0/6.0));
  double Kxx    = - 0.5*peak_number*pi*ICA*cos( peak_number * pi*(x[0] - t))/sqrt( 1.0 - ICA*sin( peak_number*pi*( x[0] - t))  );
  double traceK = Kxx/HH;
  tarch::memset(Q, .0, nVars*sizeof(double));
  Q[0]  = phi*phi*HH ;					//\tilde(\gamma)_xx
  Q[3]  = phi*phi  ;					//\tilde(\gamma)_yy
  Q[5]  = phi*phi                            ;		//\tilde(\gamma)_zz
  Q[6]  = phi*phi*(Kxx - 1.0/3.0*traceK*HH ) ;		//\tilde(A)_xx
  Q[9] =  phi*phi*(0.0 - 1.0/3.0*traceK*1.0) ;		//\tilde(A)_yy
  Q[11] = phi*phi*(0.0 - 1.0/3.0*traceK*1.0) ;		//\tilde(A)_zz
  Q[16] = sqrt(HH);				        //\alpha
  Q[13] = 2.0/(3.0*pow(HH,(5.0/3.0)))*dxH        ;	//\hat(\Gamma)^x
  Q[23] = 1.0/(2.0)*dxH*pow(HH,(-1.0/2.0))       ;		//A_x
  Q[35] = pow(HH,(-1.0/3.0))*dxH/3.0         ;		//D_xxx
  Q[38] = phi*dxphi                     ;		//D_xyy
  Q[40] = phi*dxphi                    ;		//D_xzz
  Q[53] = traceK;					//K
  Q[54] = phi;					//\phi
  Q[55] = dxphi;					//P_x
}

void applications::exahype2::ccz4::diagonal_gaugeWave(
  double * __restrict__ Q, // Q[64+0],
  const tarch::la::Vector<Dimensions,double>&  x,
  double t
) {
  constexpr int nVars = 59;
  constexpr double pi = M_PI;
  constexpr double peak_number = 2.0;
  constexpr double ICA = 0.1; ///< Amplitude of the wave, B in the note 
  double phase  = peak_number*pi*( x[0] - x[1] - t);
  double H=1+2.0*ICA*sin(phase);
  tarch::memset(Q, .0, nVars*sizeof(double));
  Q[0]  = pow(H,(-1.0/3.0))*(1+ICA*sin(phase));				//\tilde(\gamma)_xx
  Q[1]  = -pow(H,(-1.0/3.0))*(ICA*sin(phase));				//\tilde(\gamma)_xy
  Q[3]  = pow(H,(-1.0/3.0))*(1+ICA*sin(phase));				//\tilde(\gamma)_yy
  Q[5]  = pow(H,(-1.0/3.0));						//\tilde(\gamma)_zz
  Q[6]  = pow(H,(-5.0/6.0))*(peak_number*pi*ICA/2.0)*cos(phase)*(1.0 - 2.0 * (1+ICA*sin(phase)) / (3.0*H) );	//\tilde(A)_xx
  Q[7]  = pow(H,(-5.0/6.0))*(peak_number*pi*ICA/2.0)*cos(phase)*(2.0 * ICA*sin(phase) / (3.0*H) - 1.0 ); 	//\tilde(A)_xy
  Q[9]  = pow(H,(-5.0/6.0))*(peak_number*pi*ICA/2.0)*cos(phase)*(1.0 - 2.0 * (1+ICA*sin(phase)) / (3.0*H) );	//\tilde(A)_yy
  Q[11] = pow(H,(-11.0/6.0))*(peak_number*pi*ICA/3.0)*cos(phase);						//\tilde(A)_zz
  Q[13] = pow(H,(-5.0/3.0))*(4.0*peak_number*pi*ICA/3.0)*cos(phase);	//\hat(\Gamma)^x
  Q[14] = -pow(H,(-5.0/3.0))*(4.0*peak_number*pi*ICA/3.0)*cos(phase);	//\hat(\Gamma)^y
  Q[16] = 0.5*log(H);							//ln(\alpha)
  Q[23] = pow(H,-1.0)*(peak_number*pi*ICA)*cos(phase);			//A_x					
  Q[24] = -pow(H,-1.0)*(peak_number*pi*ICA)*cos(phase);			//A_y	
  Q[35] = pow(H,(-4.0/3.0))*(peak_number*pi*ICA/6.0)*cos(phase)*(1.0 + 4.0 * ICA*sin(phase) );		//D_xxx
  Q[36] = -pow(H,(-4.0/3.0))*(peak_number*pi*ICA/6.0)*cos(phase)*(3.0 + 4.0 * ICA*sin(phase) );		//D_xxy
  Q[38] = pow(H,(-4.0/3.0))*(peak_number*pi*ICA/6.0)*cos(phase)*(1.0 + 4.0 * ICA*sin(phase) );		//D_xyy
  Q[40] = -pow(H,(-4.0/3.0))*(peak_number*pi*ICA/3.0)*cos(phase);					//D_xzz
  Q[41] = -pow(H,(-4.0/3.0))*(peak_number*pi*ICA/6.0)*cos(phase)*(1.0 + 4.0 * ICA*sin(phase) );		//D_yxx
  Q[42] = pow(H,(-4.0/3.0))*(peak_number*pi*ICA/6.0)*cos(phase)*(3.0 + 4.0 * ICA*sin(phase) );		//D_yxy
  Q[41] = -pow(H,(-4.0/3.0))*(peak_number*pi*ICA/6.0)*cos(phase)*(1.0 + 4.0 * ICA*sin(phase) );		//D_yyy
  Q[46] = pow(H,(-4.0/3.0))*(peak_number*pi*ICA/3.0)*cos(phase);					//D_yzz
  Q[53] = pow(H,(-3.0/2.0))*(peak_number*pi*ICA)*cos(phase);		//K
  Q[54] = -(1.0/6.0)*log(H);						//ln(\phi)
  Q[55] = -pow(H,-1.0)*(peak_number*pi*ICA/3.0)*cos(phase);		//P_x
  Q[56] = pow(H,-1.0)*(peak_number*pi*ICA/3.0)*cos(phase);		//P_x
}

void applications::exahype2::ccz4::linearWave(
  double * __restrict__ Q, // Q[64+0],
  const tarch::la::Vector<Dimensions,double>&  X,
  double t
) {
  constexpr int nVars = 59;
  constexpr double pi = M_PI;
  constexpr double peak_number = 2.0;
  constexpr double ICA = 1e-4; ///< Amplitude of the wave, should be very small to keep linearized.
  double HH     =  ICA*sin( peak_number*pi*( X[0] - t));
  double dxHH   =  peak_number*pi*ICA*cos( peak_number * pi*(X[0] - t));
  double dtHH   = -peak_number*pi*ICA*cos( peak_number * pi*(X[0] - t));
  tarch::memset(Q, .0, nVars*sizeof(double));
  Q[0]  = 1.0  ;		//\tilde(\gamma)_xx
  Q[3]  = 1+HH ;		//\tilde(\gamma)_yy
  Q[5]  = 1-HH ;		//\tilde(\gamma)_zz
  Q[6]  = 0.0	    ;		//\tilde(A)_xx
  Q[9]  = -0.5*dtHH ;		//\tilde(A)_yy
  Q[11] = 0.5*dtHH  ;		//\tilde(A)_zz
  Q[16] = log(1.0) ;		//ln(\alpha)
  Q[35] = 0.0       ;		//D_xxx
  Q[38] = 0.5*dxHH  ;		//D_xyy
  Q[40] = -0.5*dxHH ;		//D_xzz
  Q[54] = log(1.0) ;		//ln(\phi)
}

void applications::exahype2::ccz4::flat(
  double * __restrict__ Q, // Q[64+0],
  const tarch::la::Vector<Dimensions,double>&  x,
  double t
) {
  tarch::memset(Q, .0, 59*sizeof(double));
  Q[0]  = 1.0 ;          //\tilde(\gamma)_xx
  Q[3]  = 1.0  ;          //\tilde(\gamma)_yy
  Q[5]  = 1.0                            ;    //\tilde(\gamma)_zz
  Q[16] = 1.0;
  Q[54] = 1.0;          //\phi
}

#ifdef IncludeTwoPunctures
void applications::exahype2::ccz4::ApplyTwoPunctures(
  double * __restrict__ Q, // Q[64+0],
  const tarch::la::Vector<Dimensions,double>&  X,
  double t,
  TP::TwoPunctures* tp,
  bool low_res
) {
  constexpr int nVars = 59;
  const double coor[3]={X[0],X[1],X[2]};
  double LgradQ[3*nVars];
  tarch::memset(Q, .0, nVars*sizeof(double));
  tp->Interpolate(coor,Q,low_res);				//do the interpolate
  //std::cout << coor[0] <<coor[1] << coor[2] << "\n";
  //std::cout << "real quantites without tilde" <<"\n";
  //for (int i=0;i<nVars;i++){std::cout << i <<"\t"<< Q[i] << "\n";}
  TP_bindding::SOCCZ4Cal(Q);				//calculate corresponding soccz4 quantities
  //std::cout << "after treatment" <<"\n";
  //for (int i=0;i<nVars;i++){std::cout << i <<"\t" << Q[i] << "\n";}
  TP_bindding::GradientCal(coor, Q, LgradQ, nVars, tp, low_res);	//calculate gradient for auxiliary variables
  //for (int d=0;d<3;d++)
  //for (int i=0;i<nVars;i++) {std::cout << d <<"\t" << i <<"\t" << LgradQ[d*nVars+i] << "\n";}
  TP_bindding::AuxiliaryCal(Q, LgradQ, nVars); //calculate the auxiliary variables
  //std::cout << "after treatment" <<"\n";
  //for (int i=0;i<nVars;i++){std::cout << i <<"\t" << Q[i] << "\n";}

  //exit(0);
}
#endif

