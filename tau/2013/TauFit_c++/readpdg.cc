//#define USING_NBASE31
//#define USING_NBASE39
#define USING_NBASE40
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TLine.h"
#include "TBox.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TSpline.h"
#include "TF1.h"
#include "TH2D.h"
#include "TString.h"
#include "TCut.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPaveText.h"
#include "Riostream.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
//
using namespace std;
//
enum e_coefnames {
  C_ETA3PIZUSAGE1    , // 0
  C_ETA3PIZUSAGE2    , // 1
  C_ETA3PIZUSAGE3    , // 2
  C_ETAPIMPIPPIZ     , // 3
  C_ETA3PIALL        , // 4
  C_ETANEUTRALMODES1 , // 5
  C_ETANEUTRALMODES2 , // 6
  C_ETACHARGEDMODES  , // 7
  C_KS2PIZ           , // 8
  C_KS2PIZby2        , // 9
  C_KS2PIZby2PlusHalf, //10
  C_KS2PIZSQUARED    , //11
  C_KS2PIZSQUAREPlus1, //12
  C_KS2PIby2         , //13
  C_KS2PI            , //14
  C_KS2PI_X_2PIZ_X_2 , //15
  C_OMPIMPIPPIZ      , //16
  C_OMPIMPIP         , //17
  C_OM3PIPlus2PI     , //18
  C_OM2PIZGAMMA      , //19
  C_PHI2KMKP         , //20
  C_PHI2KSKL         , //21
  C_PHI2KK_KS2PI     , //22
  C_PHI2KSKL_KS2PI   , //23
  C_PHI2KSKL_KS2PIZ  , //24
  C_Gam132_3PI3PIZ_KL, //25
  C_Gam132_ETANeuby2 , //26
  C_ONE              , //27
  C_HALF             , //28
  C_TWO              , //29
  C_ONETHIRD         , //30
  C_TWOTHIRD           //31
};
//
enum e_basegammanames {
  M_GAMMA3  , // 0
  M_GAMMA5  , // 1
  M_GAMMA9  , // 2
  M_GAMMA10 , // 3
  M_GAMMA14 , // 4
  M_GAMMA16 , // 5
  M_GAMMA20 , // 6
  M_GAMMA23 , // 7
  M_GAMMA27 , // 8
  M_GAMMA28 , // 9
  M_GAMMA30 , // 10
  M_GAMMA35 , // 11
  M_GAMMA37 , // 12
  M_GAMMA40 , // 13
  M_GAMMA42 , // 14
  M_GAMMA47 , // 15
  M_GAMMA48 , // 16
  M_GAMMA62 , // 17
  M_GAMMA70 , // 18
  M_GAMMA77 , // 19
  M_GAMMA78 , // 20
#if defined USING_NBASE31
  M_GAMMA85 , // 21
  M_GAMMA89 , // 22
#elif defined USING_NBASE39 || defined USING_NBASE40 || defined USING_NBASE47
  M_GAMMA802, // 21
  M_GAMMA803, // 22
#endif
  M_GAMMA93 , // 23
  M_GAMMA94 , // 24
  M_GAMMA104, // 25
  M_GAMMA126, // 26
  M_GAMMA128, // 27
#if defined USING_NBASE31
  M_GAMMA150, // 28
#elif defined USING_NBASE39 || defined USING_NBASE40 || defined USING_NBASE47
  M_GAMMA800, //    // 28
  M_GAMMA151, //    // 29
#endif
  M_GAMMA152, // 29 // 30
#if defined USING_NBASE39 || defined USING_NBASE40 || defined USING_NBASE47
  M_GAMMA130, //    // 31
  M_GAMMA132, //    // 32
  M_GAMMA44 , //    // 33
  M_GAMMA53 , //    // 34
  M_GAMMA49,  //    // 35
  M_GAMMA804, //    // 36
  M_GAMMA805, //    // 37
#if defined USING_NBASE40
  M_GAMMA801, //    // 38
#endif
#endif
  M_GAMMA103 //     // 39
};
//
enum e_nodegammanames {
  N_GAMMA128    ,
  N_GAMMA19BY13 ,
  N_GAMMA26BY13 ,
  N_GAMMA30     ,
  N_GAMMA76BY54 ,
  N_GAMMA152BY54,
  N_GAMMA152BY76,
  N_GAMMA16     ,
  N_GAMMA23     ,
  N_GAMMA28     ,
  N_GAMMA35     ,
  N_GAMMA40     ,
  N_GAMMA42     ,
  N_GAMMA92     ,
  N_GAMMA33     ,
  N_GAMMA106    ,
  N_GAMMA46     ,
  N_GAMMA66     ,
  N_GAMMA67     ,
  N_GAMMA20     ,
  N_GAMMA27     ,
  N_GAMMA78     ,
  N_GAMMA152    ,
  N_GAMMA76     ,
  N_GAMMA57     ,
  N_GAMMA55     ,
  N_GAMMA57BY55 ,
  N_GAMMA34     ,
  N_GAMMA39     ,
  N_GAMMA47     ,
  N_GAMMA58     ,
  N_GAMMA77     ,
  N_GAMMA8      ,
  N_GAMMA18     ,
  N_GAMMA1      ,
  N_GAMMA65     ,
  N_GAMMA75     ,
  N_GAMMA64     ,
  N_GAMMA29     ,
  N_GAMMA8BY5   ,
  N_GAMMA12     ,
  N_GAMMA25     ,
  N_GAMMA74     ,
  N_GAMMA48     ,
  N_GAMMA59     ,
  N_GAMMA60     ,
  N_GAMMA62     ,
  N_GAMMA85     ,
  N_GAMMA68     ,
  N_GAMMA69     ,
  N_GAMMA70     ,
  N_GAMMA88     ,
  N_GAMMA80     ,
  N_GAMMA80BY60 ,
  N_GAMMA81     ,
  N_GAMMA81BY69 ,
  N_GAMMA93BY60 ,
  N_GAMMA94BY69 ,
  N_GAMMA38     ,
  N_GAMMA83     ,
  N_GAMMA110    ,
  N_GAMMA89     ,
  N_GAMMA84     ,
  N_GAMMA87     ,
  N_GAMMA94     ,
  N_GAMMA3      ,
  N_GAMMA150BY66,
  N_GAMMA149    ,
  N_GAMMA5      ,
  N_GAMMA19     ,
  N_GAMMA26     ,
  N_GAMMA150    ,
  N_GAMMA2      ,
  N_GAMMA31     ,
  N_GAMMA32     ,
  N_GAMMA56     ,
  N_GAMMA63     ,
  N_GAMMA54     ,
  N_GAMMA126    ,
  N_GAMMA102    ,
  N_GAMMA79     ,
  N_GAMMA103    ,
  N_GAMMA104    ,
  N_GAMMA93     ,
  N_GAMMA82     ,
  N_GAMMA11     ,
  N_GAMMA7      ,
  N_GAMMA17     ,
  N_GAMMA37     ,
  N_GAMMA3BY5   ,
  N_GAMMA9      ,
  N_GAMMA10     ,
  N_GAMMA14     ,
  N_GAMMA13     ,
  N_GAMMA24     ,
  N_GAMMA9BY5   ,
  N_GAMMA10BY5  ,
  N_GAMMA10BY9  ,
#if defined USING_NBASE39 || defined USING_NBASE40
  N_GAMMA130    ,
  N_GAMMA132    ,
  N_GAMMA43     ,
  N_GAMMA44     ,
  N_GAMMA53     ,
  N_GAMMA800    ,
  N_GAMMA151    ,
  N_GAMMA802    ,
  N_GAMMA803    ,
  N_GAMMA136    ,
  N_GAMMA115    ,
  N_GAMMA49     ,
  N_GAMMA804    ,
  N_GAMMA805    ,
#if defined USING_NBASE40
  N_GAMMA96     ,
  N_GAMMA801    ,
#endif
#endif
  N_GAMMAALL    
};
enum P_ParNames{
  P_m_e, //0
  P_m_mu,//1
  P_tau_mu,//2
  P_m_tau,//3
  P_tau_tau,//4
  P_m_pim,//5
  P_tau_pim,//6
  P_BR_PimToMumNu,//7
  P_Delta_TauToPim_over_PimToMu,//8
  P_m_km,//9
  P_tau_km,//10
  P_BR_KmToMumNu,//11
  P_Delta_TauToKm_over_KmToMu,//12
  P_m_W,//13
  P_alpha_val,//14
  P_Delta_mu_gamma,//15
  P_Delta_tau_gamma,//16
  P_GMu,//17
  P_hbar,//18
  P_SEW,//19
  P_fpi,//20
  P_fK,//21
  P_fKfpi,//22
  P_Vud,//23
  P_Delta_Kpi,//24
  P_Delta_Rth,//25
  P_Be,//26
  P_Bmu,//27
  P_Bpi,//28
  P_BK,//29
  P_BmuBe,//30
  P_BKBpi,//31
  P_Bs,//32
  P_Ball,//33
  P_wt_Be,//34
  P_wt_Be_from_Bmu,//35
  P_wt_Be_from_tautau,//36
  P_wt_Be_from_Bpi,//37
  P_wt_Be_from_BK//38
};
const int N_ParNames = 39;
const int N_GlobalPar= 26;
//
//--- from PDG 2010 
//            MOHR     08 RVUE                2006 CODATA value
//  MOHR 2008         Reviews of Modern Physics 80 (2008) 633 
// Published in Rev.Mod.Phys.80:633-730,2008.  e-Print: arXiv:0801.0028 [physics.atom-ph]
// CODATA Recomended Values of the Fundamental Physical Constants: 2006
// http://www.slac.stanford.edu/spires/find/hep/www?eprint=arXiv:0801.0028                     
const double m_e =  0.510998910; const double e_m_e = 0.000000013 ; // MeV 
const double m_mu = 105.6583668 ; const double e_m_mu = 0.0000038 ; // MeV 
//--- from HFAG 2009
const double m_tau =  1776.7673082; const double e_m_tau = 0.1507259; // MeV 
//--- from PDG 2010 
const double tau_mu  = 2.197034e-6;const double e_tau_mu  = 0.000021e-6; // s
const double tau_tau = 290.6e-15 ; const double e_tau_tau = 1.0e-15; // s
//--- from PDG 2010 
const double m_pim = 139.57018;   const double e_m_pim   = 0.00035; // MeV
const double tau_pim = 2.6033e-8; const double e_tau_pim = 0.0005e-8; // s
const double BR_PimToMumNu = 99.98770e-2 ; const double e_BR_PimToMumNu = 0.00004e-2;
//--- from PDG 2010  [Average]
const double m_km   = 493.677; const double e_m_km = 0.013; // MeV
const double tau_km = 1.2379e-8; const double e_tau_km = 0.0021e-8; // s
const double BR_KmToMumNu = 63.60e-2; const double e_BR_KmToMumNu = 0.16e-2;
//--- from Marciano:1993sh,Decker:1994ea,Decker:1994dd
const double Delta_TauToPim_over_PimToMu = 1.0016 ; const double e_Delta_TauToPim_over_PimToMu = 0.0014; 
const double Delta_TauToKm_over_KmToMu = 1.0090 ;   const double e_Delta_TauToKm_over_KmToMu = 0.0022; 
//
// rad. corrections from to get Be from tau lifetime
// values from M. Davier, et. al., 10.1103/RevModPhys.78.1043 p.1047, arXiv:hep-ph/0507078v2 p.7, 
// Delta^L_gamma = 1 + alpha(mL)/2pi * (25/4 - pi^2)
// Delta^L_W = 1 + 3/5* m_L^2/M_W^2
//
//--- from PDG 2010
const double m_W = 80.399*1e3 ; const double e_m_W =  0.023*1e3; // MeV
//--- from PDG 2010
const double alpha_val = 7.2973525376e-3;  // at Q^2 = 0
const double alpha_err = 0.0000000050e-3;
const double Delta_mu_gamma = 1 - 42.4e-4;  const double e_Delta_mu_gamma  = alpha_err * TMath::Abs( ( 1./TMath::TwoPi()) * ( (25./4.) - TMath::Power(TMath::Pi(),2) ) );  
const double Delta_tau_gamma = 1 - 43.2e-4; const double e_Delta_tau_gamma = alpha_err * TMath::Abs( ( 1./TMath::TwoPi()) * ( (25./4.) - TMath::Power(TMath::Pi(),2) ) );  
//const double Delta_mu_W = 1 + 1.0e-6;  // recomputed later using ([m_mu,  m_W], [e_m_mu,  e_m_W])
//const double Delta_tau_W = 1 + 2.9e-4; // recomputed later using ([m_tau, m_W], [e_m_tau, e_m_W])
//
const double GMu = 1.16637e-5*1e-6 ; const double e_GMu = 0.00001e-5*1e-6; // GMu in MeV-2
const double hbar=6.58211899e-25*1e3;const double e_hbar= 0.00000016e-25*1e3; // MeV s
const double SEW = 1.0201;           const double e_SEW = 0.0003;
//
const double fpi = 132;              const double e_fpi = 2; // fpi in MeV arXiv:0706.1726 [hep-lat]
const double fK  = 157;              const double e_fK  = 2; // fK  in MeV arXiv:0706.1726 [hep-lat]
const double fKfpi= 1.189;           const double e_fKfpi=0.007;
//
const double Vud = 0.97425;          const double e_Vud = 0.00022;
//
const double Delta_Kpi=1.00034;      const double e_Delta_Kpi=0.00437;
const double Delta_Rth = 0.240;      const double e_Delta_Rth = 0.032;
//
double GlobalPar[N_GlobalPar];
double e_GlobalPar[N_GlobalPar];
// ----------------------------------------------------------------------
double MyGradientPar(TF1* func, int ipar, const int npar, const double *x, const double *par, const double *epar, double eps){
   if(eps< 1e-10 || eps > 1) {
      Warning("Derivative","parameter esp=%g out of allowed range[1e-10,1], reset to 0.01",eps);
      eps = 0.01;
   }
   //save original parameters
   double par0 = par[ipar];
   // new array
   double parv[npar];
   for (int i=0;i<npar;++i) parv[i] = par[i];
   //
   func->InitArgs(x, parv);
   double f=func->EvalPar(x,parv) ; 
   //
   double f1, f2, g1, g2, h2, d0, d2;
   double h = eps*epar[ipar];
   parv[ipar] = par0 + h;     f1 = func->EvalPar(x,parv);
   parv[ipar] = par0 - h;     f2 = func->EvalPar(x,parv);
   parv[ipar] = par0 + h/2;   g1 = func->EvalPar(x,parv);
   parv[ipar] = par0 - h/2;   g2 = func->EvalPar(x,parv);
   //compute the central differences
   h2    = 1/(2.*h);
   d0    = f1 - f2;
   d2    = 2*(g1 - g2);
   double  grad = h2*(4*d2 - d0)/3.;
//   //
//   cout << "eps = " << eps  << " epar = " << epar[ipar] << " h = " << h << endl;   
//   cout << "par0     = " << par0   << " f  = " << f  << endl;
//   cout << "par0+h   = " << par0+h << " f1 = " << f1 << endl;
//   cout << "par0-h   = " << par0-h << " f2 = " << f2 << endl;
//   cout << "par0+h/2 = " << par0+h << " g1 = " << g1 << endl;
//   cout << "par0-h/2 = " << par0-h << " g2 = " << g2 << endl;
//   //
   return grad;
}
// ----------------------------------------------------------------------
void calc_generic(double *par, double *epar, double** CovarianceMatrix, double (*fcn)(double*, double*), double& func1_val, double& func1_err, double* func1_evec){ 
  //
  int ipar,jpar;
  //
  const double eps=1e-2;
  //
  const double xdummy[1]={0};
  TF1* func1 = new TF1 ("func1", fcn, xdummy[0], xdummy[0]+1, N_ParNames);
  func1->InitArgs(xdummy,par);
  func1_val = func1->EvalPar(xdummy,par);
  //
  double* deriv = new double[N_ParNames];
  for (ipar=0;ipar<N_ParNames;++ipar) {
    deriv[ipar] = (epar[ipar]) ? MyGradientPar(func1,ipar,N_ParNames,xdummy,par,epar,eps) : 0;
  }
  //
  double func1_err2=0;
  for (ipar=0;ipar<N_ParNames;++ipar) {
    func1_evec[ipar] = deriv[ipar] * epar[ipar]; 
    for (jpar=0;jpar<N_ParNames;++jpar) {
      func1_err2 += deriv[ipar] * CovarianceMatrix[ipar][jpar] * deriv[jpar]; 
    }
  }
  func1_err = TMath::Sqrt(func1_err2);
  //
  delete [] deriv;
  delete func1;
}
// ----------------------------------------------------------------------
double func_Vus_Kpi(double *x, double * par){
  double dummy = x[0];
  double value = TMath::Sqrt( 
			     par[P_BKBpi] * par[P_Vud]*par[P_Vud] * (1./(par[P_fKfpi]*par[P_fKfpi])) *
			     (TMath::Power(1 - TMath::Power(par[P_m_pim]/par[P_m_tau],2),2) / TMath::Power(1 - TMath::Power(par[P_m_km]/par[P_m_tau],2),2)) *
			     (1./par[P_Delta_Kpi])); 
  return value;
}
// ----------------------------------------------------------------------
double func_Vus_K(double *x, double * par){
  double dummy = x[0];
  double value = TMath::Sqrt( 
			     par[P_BK] * 16.0 * TMath::Pi() * par[P_hbar] * 1./ 
			     (TMath::Power(par[P_GMu],2) * TMath::Power(par[P_fK],2) * TMath::Power(par[P_m_tau],3) * par[P_tau_tau] * 
			      TMath::Power(1 - TMath::Power(par[P_m_km]/par[P_m_tau],2),2) * par[P_SEW]) );
  return value;
}
// ----------------------------------------------------------------------
double func_fmufe(double *x, double * par){
  double dummy = x[0];
  double ratio1 = (par[P_m_e]*par[P_m_e]) / (par[P_m_tau]*par[P_m_tau]);
  double value1 = 1 - (8 * ratio1) + (8 * ratio1 * ratio1 * ratio1) - (ratio1 * ratio1 * ratio1 * ratio1) - (12 * ratio1 * ratio1 * log(ratio1));
  double ratio2 = (par[P_m_mu]*par[P_m_mu]) / (par[P_m_tau]*par[P_m_tau]);
  double value2 = 1 - (8 * ratio2) + (8 * ratio2 * ratio2 * ratio2) - (ratio2 * ratio2 * ratio2 * ratio2) - (12 * ratio2 * ratio2 * log(ratio2));
  double fmufe = value2/value1;
  return fmufe;
}
// ----------------------------------------------------------------------
double func_gmuge(double *x, double * par){
  double dummy = x[0];
  double fmufe = func_fmufe(x,par);
  double gmuge = TMath::Sqrt(par[P_BmuBe]/fmufe);
  return gmuge;
}
// ----------------------------------------------------------------------
double func_Be_from_Bmu(double *x, double * par){
  double dummy = x[0];
  double fmufe = func_fmufe(x,par);
  double Be_from_Bmu = par[P_Bmu]/fmufe;
  return Be_from_Bmu;
}
// ----------------------------------------------------------------------
double func_Be_from_Bpi(double *x, double * par){
  double dummy = x[0];
/*
(gtau/gmu)^2 = (tau_mu/tau_tau) * (m_mu/m_tau)^5 * Be * [f(m_e^2/m_mu^2) / f(m_e^2/m_tau^2)] * (Delta_mu_W/Delta_tau_W) * (Delta_mu_gamma/Delta_tau_gamma)
(gtau/gmu)^2 = (BR_TauToPiNu / BR_PimToMumNu) * (tau_pi/tau_tau) * (2*m_pi*m_mu^2/m_tau^3) * (1/.pow((1-pow(M_Pi/M_Tau,2))/(1-pow(M_Mu/M_Pi,2)),2)) *  (1/Delta_Pi)
=>
 Be = (BR_TauToPiNu / BR_PimToMumNu) * (tau_pi/tau_mu) * (2*m_pi *m_tau^2/m_mu^3) * (1/.pow((1-pow(m_pim/m_tau,2))/(1-pow(m_mu/m_pim,2)),2)) *  (1/Delta_TauToPim_over_PimToMu)
      * 1/[f(m_e^2/m_mu^2) / f(m_e^2/m_tau^2)] * 1/[(Delta_mu_W/Delta_tau_W)] * 1/[(Delta_mu_gamma/Delta_tau_gamma)]
*/
  double ratio1 = (par[P_m_e]*par[P_m_e]) / (par[P_m_tau]*par[P_m_tau]);
  double value_etau = 1 - (8 * ratio1) + (8 * ratio1 * ratio1 * ratio1) - (ratio1 * ratio1 * ratio1 * ratio1) - (12 * ratio1 * ratio1 * log(ratio1));
  double ratio2 = (par[P_m_e]*par[P_m_e]) / (par[P_m_mu]*par[P_m_mu]);
  double value_emu  = 1 - (8 * ratio2) + (8 * ratio2 * ratio2 * ratio2) - (ratio2 * ratio2 * ratio2 * ratio2) - (12 * ratio2 * ratio2 * log(ratio2));
  // Delta^L_W = 1 + 3/5* m_L^2/M_W^2
  double Delta_tau = 1 + 0.6 * (par[P_m_tau]*par[P_m_tau]) / (par[P_m_W]*par[P_m_W]);
  double Delta_mu  = 1 + 0.6 * (par[P_m_mu]*par[P_m_mu]) / (par[P_m_W]*par[P_m_W]);
  return 
    ( par[P_Bpi] / par[P_BR_PimToMumNu] ) * 
    ( par[P_tau_pim] / par[P_tau_mu] ) *
    (( 2 * par[P_m_pim] * TMath::Power(par[P_m_tau],2)) / TMath::Power(par[P_m_mu],3)) *
    TMath::Power((1-TMath::Power(par[P_m_mu]/par[P_m_pim],2))/(1-TMath::Power(par[P_m_pim]/par[P_m_tau],2)),2) *
    (1/par[P_Delta_TauToPim_over_PimToMu]) *
    (value_etau/value_emu) *
    (Delta_tau/Delta_mu) *
    (par[P_Delta_tau_gamma]/par[P_Delta_mu_gamma]);
}
// ----------------------------------------------------------------------
double func_Be_from_BK(double *x, double * par){
  double dummy = x[0];
/*
(gtau/gmu)^2 = (tau_mu/tau_tau) * (m_mu/m_tau)^5 * Be * [f(m_e^2/m_mu^2) / f(m_e^2/m_tau^2)] * (Delta_mu_W/Delta_tau_W) * (Delta_mu_gamma/Delta_tau_gamma)
(gtau/gmu)^2 = (BR_TauToKNu / BR_KmToMumNu) * (tau_K/tau_tau) * (2*m_K*m_mu^2/m_tau^3) * (1/.pow((1-pow(M_K/M_Tau,2))/(1-pow(M_Mu/M_K,2)),2)) *  (1/Delta_K)
=>
 Be = (BR_TauToKNu / BR_KmToMumNu) * (tau_K/tau_mu) * (2*m_K *m_tau^2/m_mu^3) * (1/.pow((1-pow(m_km/m_tau,2))/(1-pow(m_mu/m_km,2)),2)) *  (1/Delta_TauToKm_over_KmToMu)
      * 1/[f(m_e^2/m_mu^2) / f(m_e^2/m_tau^2)] * 1/[(Delta_mu_W/Delta_tau_W)] * 1/[(Delta_mu_gamma/Delta_tau_gamma)]
*/
  double ratio1 = (par[P_m_e]*par[P_m_e]) / (par[P_m_tau]*par[P_m_tau]);
  double value_etau = 1 - (8 * ratio1) + (8 * ratio1 * ratio1 * ratio1) - (ratio1 * ratio1 * ratio1 * ratio1) - (12 * ratio1 * ratio1 * log(ratio1));
  double ratio2 = (par[P_m_e]*par[P_m_e]) / (par[P_m_mu]*par[P_m_mu]);
  double value_emu  = 1 - (8 * ratio2) + (8 * ratio2 * ratio2 * ratio2) - (ratio2 * ratio2 * ratio2 * ratio2) - (12 * ratio2 * ratio2 * log(ratio2));
  // Delta^L_W = 1 + 3/5* m_L^2/M_W^2
  double Delta_tau = 1 + 0.6 * (par[P_m_tau]*par[P_m_tau]) / (par[P_m_W]*par[P_m_W]);
  double Delta_mu  = 1 + 0.6 * (par[P_m_mu]*par[P_m_mu]) / (par[P_m_W]*par[P_m_W]);
  return 
    ( par[P_BK] / par[P_BR_KmToMumNu] ) * 
    ( par[P_tau_km] / par[P_tau_mu] ) *
    (( 2 * par[P_m_km] * TMath::Power(par[P_m_tau],2)) / TMath::Power(par[P_m_mu],3)) *
    TMath::Power((1-TMath::Power(par[P_m_mu]/par[P_m_km],2))/(1-TMath::Power(par[P_m_km]/par[P_m_tau],2)),2) *
    (1/par[P_Delta_TauToKm_over_KmToMu]) *
    (value_etau/value_emu) *
    (Delta_tau/Delta_mu) *
    (par[P_Delta_tau_gamma]/par[P_Delta_mu_gamma]);
}
// ----------------------------------------------------------------------
double func_Bpi_univ(double *x, double * par){
  double dummy = x[0];
  // E. Gamiz, et. al., http://arxiv.org/pdf/0709.0282v1 Eqn 4.1
  // BR_TauToPiNu = BR_PimToMumNu * pow(M_Tau,3) * Tau_Tau * 1./(2 * pow(M_Mu,2) * M_Pi * Tau_Pi) * pow((1-pow(M_Pi/M_Tau,2))/(1-pow(M_Mu/M_Pi,2)),2) * Delta_Pi ;
  return 
    par[P_BR_PimToMumNu] * 
    TMath::Power(par[P_m_tau],3) * 
    par[P_tau_tau] * 
    1/(2 * TMath::Power(par[P_m_mu],2) * par[P_m_pim] * par[P_tau_pim]) * 
    TMath::Power((1-TMath::Power(par[P_m_pim]/par[P_m_tau],2))/(1-TMath::Power(par[P_m_mu]/par[P_m_pim],2)),2) * 
    par[P_Delta_TauToPim_over_PimToMu];
}
// ----------------------------------------------------------------------
double func_gtaugmu_pi(double *x, double * par){
  double dummy = x[0];
  double Bpi_univ = func_Bpi_univ(x,par);
  double gtaugmu_pi = sqrt(par[P_Bpi]/Bpi_univ);
  return gtaugmu_pi;
}
// ----------------------------------------------------------------------
double func_BK_univ(double *x, double * par){
  double dummy = x[0];
  // E. Gamiz, et. al., http://arxiv.org/pdf/0709.0282v1 Eqn 4.1
  // BR_TauToPiNu = BR_PimToMumNu * pow(M_Tau,3) * Tau_Tau * 1./(2 * pow(M_Mu,2) * M_Pi * Tau_Pi) * pow((1-pow(M_Pi/M_Tau,2))/(1-pow(M_Mu/M_Pi,2)),2) * Delta_Pi ;
  return 
    par[P_BR_KmToMumNu] * 
    TMath::Power(par[P_m_tau],3) * 
    par[P_tau_tau] * 
    1/(2 * TMath::Power(par[P_m_mu],2) * par[P_m_km] * par[P_tau_km]) * 
    TMath::Power((1-TMath::Power(par[P_m_km]/par[P_m_tau],2))/(1-TMath::Power(par[P_m_mu]/par[P_m_km],2)),2) * 
    par[P_Delta_TauToKm_over_KmToMu];
}
// ----------------------------------------------------------------------
double func_gtaugmu_K(double *x, double * par){
  double dummy = x[0];
  double BK_univ = func_BK_univ(x,par);
  double gtaugmu_K = sqrt(par[P_BK]/BK_univ);
  return gtaugmu_K;
}
// ----------------------------------------------------------------------
double func_Be_from_tautau(double *x, double * par){
  // double Be_from_tautau = (phspf_mebymtau/phspf_mebymmu) * TMath::Power((m_tau/m_mu),5) * (tau_tau/tau_mu) * (Delta_tau_W/Delta_mu_W) * (Delta_tau_gamma/Delta_mu_gamma);
  double dummy = x[0];
  double ratio1 = (par[P_m_e]*par[P_m_e]) / (par[P_m_tau]*par[P_m_tau]);
  double value1 = 1 - (8 * ratio1) + (8 * ratio1 * ratio1 * ratio1) - (ratio1 * ratio1 * ratio1 * ratio1) - (12 * ratio1 * ratio1 * log(ratio1));
  double ratio2 = (par[P_m_e]*par[P_m_e]) / (par[P_m_mu]*par[P_m_mu]);
  double value2 = 1 - (8 * ratio2) + (8 * ratio2 * ratio2 * ratio2) - (ratio2 * ratio2 * ratio2 * ratio2) - (12 * ratio2 * ratio2 * log(ratio2));
  // Delta^L_W = 1 + 3/5* m_L^2/M_W^2
  double Delta_mu_gamma = 1 + 0.6 * (par[P_m_mu]*par[P_m_mu]) / (par[P_m_W] * par[P_m_W]);
  double Delta_tau_gamma = 1 + 0.6 * (par[P_m_tau]*par[P_m_tau]) / (par[P_m_W] * par[P_m_W]);
  return 
    (value1/value2) * 
    TMath::Power(par[P_m_tau]/par[P_m_mu],5) * 
    (par[P_tau_tau]/par[P_tau_mu]) * 
    (Delta_tau_gamma/Delta_mu_gamma) * 
    (par[P_Delta_tau_gamma]/par[P_Delta_mu_gamma]) ;
}
// ----------------------------------------------------------------------
double func_gtaugmu_e(double *x, double * par){
  double dummy = x[0];
  double Be_from_tautau = func_Be_from_tautau(x,par);
  double gtaugmu_e = sqrt(par[P_Be]/Be_from_tautau);
  return gtaugmu_e;
}
// ----------------------------------------------------------------------
double func_Bmu_from_tautau(double *x, double * par){
  // double Bmu_from_tautau = (phspf_mmubymtau/phspf_mebymmu) * TMath::Power((m_tau/m_mu),5) * (tau_tau/tau_mu) * (Delta_tau_W/Delta_mu_W) * (Delta_tau_gamma/Delta_mu_gamma);
  double dummy = x[0];
  double ratio1 = (par[P_m_mu]*par[P_m_mu]) / (par[P_m_tau]*par[P_m_tau]);
  double value1 = 1 - (8 * ratio1) + (8 * ratio1 * ratio1 * ratio1) - (ratio1 * ratio1 * ratio1 * ratio1) - (12 * ratio1 * ratio1 * log(ratio1));
  double ratio2 = (par[P_m_e]*par[P_m_e]) / (par[P_m_mu]*par[P_m_mu]);
  double value2 = 1 - (8 * ratio2) + (8 * ratio2 * ratio2 * ratio2) - (ratio2 * ratio2 * ratio2 * ratio2) - (12 * ratio2 * ratio2 * log(ratio2));
  // Delta^L_W = 1 + 3/5* m_L^2/M_W^2
  double Delta_mu_gamma = 1 + 0.6 * (par[P_m_mu]*par[P_m_mu]) / (par[P_m_W] * par[P_m_W]);
  double Delta_tau_gamma = 1 + 0.6 * (par[P_m_tau]*par[P_m_tau]) / (par[P_m_W] * par[P_m_W]);
  return 
    (value1/value2) * 
    TMath::Power(par[P_m_tau]/par[P_m_mu],5) * 
    (par[P_tau_tau]/par[P_tau_mu]) * 
    (Delta_tau_gamma/Delta_mu_gamma) * 
    (par[P_Delta_tau_gamma]/par[P_Delta_mu_gamma]) ;
}
// ----------------------------------------------------------------------
double func_gtauge_mu(double *x, double * par){
  double dummy = x[0];
  double Bmu_from_tautau = func_Bmu_from_tautau(x,par);
  double gtauge_mu = sqrt(par[P_Bmu]/Bmu_from_tautau);
  return gtauge_mu;
}
// ----------------------------------------------------------------------
void calc_average(const int n, const double* value, double** covmat, double& average, double& error, double* wt){
  //
  int i,j;
  //
  TMatrixD covmatrix(n,n);
  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      covmatrix[i][j] = covmat[i][j];
    }
  }
  //
  TMatrixD inverse = covmatrix;
  double determinant;
  inverse.Invert(&determinant);
  //
  for (i=0;i<n;++i) {
    wt[i] = 0;
    for (j=0;j<n;++j) {
      wt[i] += inverse[i][j];
    }
  }
  //
  double sumwt=0;
  for (i=0;i<n;++i) {
    sumwt += wt[i];
  }
  //
  for (i=0;i<n;++i) {
    wt[i] /= sumwt;
  }
  //
  average = 0;
  for (i=0;i<n;++i) {
    average += wt[i] * value[i];
  }
  error = sqrt(1/sumwt);
  //
  double pull[n];
  for (i=0;i<n;++i) {
    pull[i] = value[i] - average;
  }
  //
  double chi2 = 0;
  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      chi2 += pull[i] * inverse[i][j] * pull[j];
    }
  }
  //  cout << "chi2 = " << chi2 << " n = " << n << " Prob = " << TMath::Prob(chi2,n-1) << endl;
}
// ----------------------------------------------------------------------
double func_Bhadrons(double* x, double* par) {
  //
  double dummy=x[0];
  //
  double Be = par[P_Be];
  double Be_from_Bmu = func_Be_from_Bmu(x,par);
  double Be_from_tautau = func_Be_from_tautau(x,par);
  double Be_from_Bpi = func_Be_from_Bpi(x,par);
  double Be_from_BK = func_Be_from_BK(x,par);
  double Btotal = par[P_Ball];
  //
  //  cout << "Be = " << Be << " Be_from_Bmu = " << Be_from_Bmu << " Be_from_tautau = " << Be_from_tautau << " Be_from_Bpi = " << Be_from_Bpi << " Be_from_BK  = " << Be_from_BK  << endl;
  double Be_univ_wt[5] = {par[P_wt_Be], par[P_wt_Be_from_Bmu], par[P_wt_Be_from_tautau], par[P_wt_Be_from_Bpi], par[P_wt_Be_from_BK]};
  //cout << Be_univ_wt[0] << " " << Be_univ_wt[0]  << " " << Be_univ_wt[2] << " " << Be_univ_wt[3] << " " << Be_univ_wt[4]  << endl;
  double Be_univ = Be_univ_wt[0] * Be + Be_univ_wt[1] * Be_from_Bmu + Be_univ_wt[2] * Be_from_tautau + Be_univ_wt[3] * Be_from_Bpi + Be_univ_wt[4] * Be_from_BK;
  //
  double fmufe = func_fmufe(x,par);
  double Bhadrons = Btotal - (1. + fmufe) * Be_univ;
  //cout << "fmufe = " << fmufe << " Be_univ = " << Be_univ << " Bhadrons = " << Bhadrons << endl;
  return Bhadrons;
}
// ----------------------------------------------------------------------
double func_Rhadrons(double* x, double* par){
  //
  double dummy=x[0];
  //
  double Be = par[P_Be];
  double Be_from_Bmu = func_Be_from_Bmu(x,par);
  double Be_from_tautau = func_Be_from_tautau(x,par);
  double Be_from_Bpi = func_Be_from_Bpi(x,par);
  double Be_from_BK = func_Be_from_BK(x,par);
  double Btotal = par[P_Ball];
  //
  double Be_univ_wt[5] = {par[P_wt_Be], par[P_wt_Be_from_Bmu], par[P_wt_Be_from_tautau], par[P_wt_Be_from_Bpi], par[P_wt_Be_from_BK]};
  double Be_univ = Be_univ_wt[0] * Be + Be_univ_wt[1] * Be_from_Bmu + Be_univ_wt[2] * Be_from_tautau + Be_univ_wt[3] * Be_from_Bpi + Be_univ_wt[4] * Be_from_BK;
  //
  double fmufe = func_fmufe(x,par);
  double Bhadrons = Btotal - (1. + fmufe) * Be_univ;
  double Rhadrons = Bhadrons/Be_univ;
  return Rhadrons;
}
// ----------------------------------------------------------------------
double func_Rstrange(double* x, double* par){
  double dummy=x[0];
  //
  double Be = par[P_Be];
  double Be_from_Bmu = func_Be_from_Bmu(x,par);
  double Be_from_tautau = func_Be_from_tautau(x,par);
  double Be_from_Bpi = func_Be_from_Bpi(x,par);
  double Be_from_BK = func_Be_from_BK(x,par);
  //
  double Be_univ_wt[5] = {par[P_wt_Be], par[P_wt_Be_from_Bmu], par[P_wt_Be_from_tautau], par[P_wt_Be_from_Bpi], par[P_wt_Be_from_BK]};
  double Be_univ = Be_univ_wt[0] * Be + Be_univ_wt[1] * Be_from_Bmu + Be_univ_wt[2] * Be_from_tautau + Be_univ_wt[3] * Be_from_Bpi + Be_univ_wt[4] * Be_from_BK;
  //
  double Bstrange = par[P_Bs];
  double Rstrange = Bstrange/Be_univ;
  return Rstrange;
}
// ----------------------------------------------------------------------
double func_Rnonstrange(double* x, double* par){
  double dummy=x[0];
  //
  double Be = par[P_Be];
  double Be_from_Bmu = func_Be_from_Bmu(x,par);
  double Be_from_tautau = func_Be_from_tautau(x,par);
  double Be_from_Bpi = func_Be_from_Bpi(x,par);
  double Be_from_BK = func_Be_from_BK(x,par);
  double Btotal = par[P_Ball];
  //
  double Be_univ_wt[5] = {par[P_wt_Be], par[P_wt_Be_from_Bmu], par[P_wt_Be_from_tautau], par[P_wt_Be_from_Bpi], par[P_wt_Be_from_BK]};
  double Be_univ = Be_univ_wt[0] * Be + Be_univ_wt[1] * Be_from_Bmu + Be_univ_wt[2] * Be_from_tautau + Be_univ_wt[3] * Be_from_Bpi + Be_univ_wt[4] * Be_from_BK;
  //
  double fmufe = func_fmufe(x,par);
  double Bhadrons = Btotal - (1. + fmufe) * Be_univ;
  double Rhadrons = Bhadrons/Be_univ;
  //
  double Bstrange = par[P_Bs];
  double Rstrange = Bstrange/Be_univ;
  //
  double Rnonstrange = Rhadrons - Rstrange;
  return Rnonstrange;
}
// ----------------------------------------------------------------------
double func_Vus_strange(double* x, double* par){
  double dummy=x[0];
  //
  double Be = par[P_Be];
  double Bmu = par[P_Bmu];
  double Btotal = par[P_Ball];
  //
  double fmufe = func_fmufe(x,par);
  double Be_from_Bmu = Bmu/fmufe;
  double Be_from_tautau = func_Be_from_tautau(x,par);
  double Be_from_Bpi = func_Be_from_Bpi(x,par);
  double Be_from_BK = func_Be_from_BK(x,par);
  //
  double Be_univ_wt[5] = {par[P_wt_Be], par[P_wt_Be_from_Bmu], par[P_wt_Be_from_tautau], par[P_wt_Be_from_Bpi], par[P_wt_Be_from_BK]};
  double Be_univ = Be_univ_wt[0] * Be + Be_univ_wt[1] * Be_from_Bmu + Be_univ_wt[2] * Be_from_tautau + Be_univ_wt[3] * Be_from_Bpi + Be_univ_wt[4] * Be_from_BK;
  //
  double Bhadrons = Btotal - (1. + fmufe) * Be_univ; 
  double Rhadrons = Bhadrons/Be_univ;
  //
  double Bstrange = par[P_Bs];
  double Rstrange = Bstrange/Be_univ;
  //
  double Rnonstrange = Rhadrons - Rstrange;
  //
  double Vud = par[P_Vud];
  double Delta_Rth = par[P_Delta_Rth];
  //
  double Vus_strange = TMath::Sqrt(Rstrange / ( Rnonstrange/(Vud*Vud) - Delta_Rth));
  return Vus_strange;
}
// ----------------------------------------------------------------------
void calc_results(FILE* thisfile, 
		  string* basetitle, double* basevalue_fit, double* baseerror_fit, double** basecov_fit, double** basecorr_fit, 
		  double* NodeValue, double* NodeError, double** NodeErrorMatrix) {
  int i;
  //
  // Correlation between G(tau- --> e- nubar(e) nu(tau)) and G(tau- --> mu- nubar(mu) nu(tau))
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA5].data(),basetitle[M_GAMMA3].data(),
	  100.*NodeErrorMatrix[N_GAMMA5][N_GAMMA3]/sqrt(NodeErrorMatrix[N_GAMMA5][N_GAMMA5]*NodeErrorMatrix[N_GAMMA3][N_GAMMA3]));
  //
  // Correlation between G(tau- --> e- nubar(e) nu(tau)) and G(tau- --> pi- nu(tau))
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA5].data(),basetitle[M_GAMMA9].data(),
	  100.*NodeErrorMatrix[N_GAMMA5][N_GAMMA9]/sqrt(NodeErrorMatrix[N_GAMMA5][N_GAMMA5]*NodeErrorMatrix[N_GAMMA9][N_GAMMA9]));
  //
  // Correlation between G(tau- --> e- nubar(e) nu(tau)) and G(tau- --> K- nu(tau))
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA5].data(),basetitle[M_GAMMA10].data(),
	  100.*NodeErrorMatrix[N_GAMMA5][N_GAMMA10]/sqrt(NodeErrorMatrix[N_GAMMA5][N_GAMMA5]*NodeErrorMatrix[N_GAMMA10][N_GAMMA10]));
  //
  // Correlation between G(tau- --> e- nubar(e) nu(tau)) and G(tau- --> strange)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA5].data(),"strange nu(tau)",
	  100.*NodeErrorMatrix[N_GAMMA5][N_GAMMA110]/sqrt(NodeErrorMatrix[N_GAMMA5][N_GAMMA5]*NodeErrorMatrix[N_GAMMA110][N_GAMMA110]));
  //
  // Correlation between G(tau- --> e- nubar(e) nu(tau)) and G(tau- --> all)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA5].data(),"sum_of_all_decays",
	  (NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]==0) ? 0 : 
	  100.*NodeErrorMatrix[N_GAMMA5][N_GAMMAALL]/sqrt(NodeErrorMatrix[N_GAMMA5][N_GAMMA5]*NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]));
  //
  // Correlation between G(tau- --> mu- nubar(mu) nu(tau)) and G(tau- --> pi- nu(tau))
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA3].data(),basetitle[M_GAMMA9].data(),
	  100.*NodeErrorMatrix[N_GAMMA3][N_GAMMA9]/sqrt(NodeErrorMatrix[N_GAMMA3][N_GAMMA3]*NodeErrorMatrix[N_GAMMA9][N_GAMMA9]));
  //
  // Correlation between G(tau- --> mu- nubar(mu) nu(tau)) and G(tau- --> K- nu(tau))
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA3].data(),basetitle[M_GAMMA10].data(),
	  100.*NodeErrorMatrix[N_GAMMA3][N_GAMMA10]/sqrt(NodeErrorMatrix[N_GAMMA3][N_GAMMA3]*NodeErrorMatrix[N_GAMMA10][N_GAMMA10]));
  //
  // Correlation between G(tau- --> mu- nubar(mu) nu(tau)) and G(tau- --> strange)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA3].data(),"strange nu(tau)",
	  100.*NodeErrorMatrix[N_GAMMA3][N_GAMMA110]/sqrt(NodeErrorMatrix[N_GAMMA3][N_GAMMA3]*NodeErrorMatrix[N_GAMMA110][N_GAMMA110]));
  //
  // Correlation between G(tau- --> mu- nubar(mu) nu(tau)) and G(tau- --> all)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA3].data(),"sum_of_all_decays",
	  (NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]==0) ? 0 : 
	  100.*NodeErrorMatrix[N_GAMMA3][N_GAMMAALL]/sqrt(NodeErrorMatrix[N_GAMMA3][N_GAMMA3]*NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]));
  //
  // Correlation between G(tau- --> pi- nu(tau)) and G(tau- --> K- nu(tau))
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA9].data(),basetitle[M_GAMMA10].data(),
	  100.*NodeErrorMatrix[N_GAMMA9][N_GAMMA10]/sqrt(NodeErrorMatrix[N_GAMMA9][N_GAMMA9]*NodeErrorMatrix[N_GAMMA10][N_GAMMA10]));
  //
  // Correlation between G(tau- --> pi- nu(tau)) and G(tau- --> strange)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA9].data(),"strange nu(tau)",
	  100.*NodeErrorMatrix[N_GAMMA9][N_GAMMA110]/sqrt(NodeErrorMatrix[N_GAMMA9][N_GAMMA9]*NodeErrorMatrix[N_GAMMA110][N_GAMMA110]));
  //
  // Correlation between G(tau- --> pi- nu(tau)) and G(tau- --> all)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA9].data(),"sum_of_all_decays",
	  (NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]==0) ? 0 : 
	  100.*NodeErrorMatrix[N_GAMMA9][N_GAMMAALL]/sqrt(NodeErrorMatrix[N_GAMMA9][N_GAMMA9]*NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]));
  //
  // Correlation between G(tau- --> K- nu(tau)) and G(tau- --> strange)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA10].data(),"strange nu(tau)",
	  100.*NodeErrorMatrix[N_GAMMA10][N_GAMMA110]/sqrt(NodeErrorMatrix[N_GAMMA10][N_GAMMA10]*NodeErrorMatrix[N_GAMMA110][N_GAMMA110]));
  //
  // Correlation between G(tau- --> K- nu(tau)) and G(tau- --> all)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA10].data(),"sum_of_all_decays",
	  (NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]==0) ? 0 : 
	  100.*NodeErrorMatrix[N_GAMMA10][N_GAMMAALL]/sqrt(NodeErrorMatrix[N_GAMMA10][N_GAMMA10]*NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]));
  //
  // Correlation between G(tau- --> strange) and G(tau- --> all)
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  "strange nu(tau)","sum_of_all_decays",
	  (NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]==0) ? 0 : 
	  100.*NodeErrorMatrix[N_GAMMA110][N_GAMMAALL]/sqrt(NodeErrorMatrix[N_GAMMA110][N_GAMMA110]*NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]));
  //
  // Correlation between G(tau- --> eta K- pi0 nu(tau)) and G(tau- --> eta pi- Kbar0 nu(tau))
  //
#if defined USING_NBASE39 || defined USING_NBASE40
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n",
	  basetitle[M_GAMMA130].data(),basetitle[M_GAMMA132].data(),
	  100.*NodeErrorMatrix[N_GAMMA130][N_GAMMA132]/sqrt(NodeErrorMatrix[N_GAMMA130][N_GAMMA130]*NodeErrorMatrix[N_GAMMA132][N_GAMMA132]));
#endif
  //
  const double   Be = basevalue_fit[M_GAMMA5];
  const double e_Be = baseerror_fit[M_GAMMA5];
  const double   Bmu = basevalue_fit[M_GAMMA3];
  const double e_Bmu = baseerror_fit[M_GAMMA3];
  const double   BmuBe = NodeValue[N_GAMMA3BY5];
  const double e_BmuBe = NodeError[N_GAMMA3BY5];
  const double   Bpi = basevalue_fit[M_GAMMA9];
  const double e_Bpi = baseerror_fit[M_GAMMA9];
  const double   BK  = basevalue_fit[M_GAMMA10];
  const double e_BK  = baseerror_fit[M_GAMMA10];
  const double   BKBpi = NodeValue[N_GAMMA10BY9];
  const double e_BKBpi = NodeError[N_GAMMA10BY9];
  double Bstrange = NodeValue[N_GAMMA110];
  double e_Bstrange = NodeError[N_GAMMA110];
  double Btotal = NodeValue[N_GAMMAALL];
  double e_Btotal = NodeError[N_GAMMAALL];
  //
  int ipar,jpar;
  //
  double *par  = new double [N_ParNames];
  double *epar = new double [N_ParNames];
  //
  for (ipar=0;ipar<N_ParNames;++ipar) {
    if (ipar<N_GlobalPar) { 
      par[ipar] = GlobalPar[ipar]; 
      epar[ipar] = e_GlobalPar[ipar];
    } else {
      par[ipar] = 0;
      epar[ipar] = 0;
    }
  }
  //
  par[P_Be]    = NodeValue[N_GAMMA5];     epar[P_Be]    = sqrt(NodeErrorMatrix[N_GAMMA5][N_GAMMA5]);
  par[P_Bmu]   = NodeValue[N_GAMMA3];     epar[P_Bmu]   = sqrt(NodeErrorMatrix[N_GAMMA3][N_GAMMA3]);
  par[P_Bpi]   = NodeValue[N_GAMMA9];     epar[P_Bpi]   = sqrt(NodeErrorMatrix[N_GAMMA9][N_GAMMA9]);
  par[P_BK]    = NodeValue[N_GAMMA10];    epar[P_BK]    = sqrt(NodeErrorMatrix[N_GAMMA10][N_GAMMA10]);
  par[P_BmuBe] = NodeValue[N_GAMMA3BY5];  epar[P_BmuBe] = sqrt(NodeErrorMatrix[N_GAMMA3BY5][N_GAMMA3BY5]);
  par[P_BKBpi] = NodeValue[N_GAMMA10BY9]; epar[P_BKBpi] = sqrt(NodeErrorMatrix[N_GAMMA10BY9][N_GAMMA10BY9]);
  par[P_Bs]    = NodeValue[N_GAMMA110];   epar[P_Bs]    = sqrt(NodeErrorMatrix[N_GAMMA110][N_GAMMA110]);
  par[P_Ball]  = NodeValue[N_GAMMAALL];   epar[P_Ball]  = sqrt(NodeErrorMatrix[N_GAMMAALL][N_GAMMAALL]);
  //  cout << par[P_Ball] << " " <<  epar[P_Ball]  << endl;
  //
  double** CovarianceMatrix = new double*[N_ParNames]; for (ipar=0;ipar<N_ParNames;++ipar) CovarianceMatrix[ipar] = new double[N_ParNames];
  //
  for (ipar=0;ipar<N_ParNames;++ipar) {
    for (jpar=0;jpar<N_ParNames;++jpar) {
      if (ipar==jpar) {
	CovarianceMatrix[ipar][ipar] = epar[ipar]*epar[ipar];
      } else {
	CovarianceMatrix[ipar][jpar] = 0;
      }
    }
  }
  //
  CovarianceMatrix[P_Delta_mu_gamma][P_Delta_tau_gamma] = CovarianceMatrix[P_Delta_tau_gamma][P_Delta_mu_gamma] = epar[P_Delta_mu_gamma]*epar[P_Delta_tau_gamma];
  //
  CovarianceMatrix[P_Be][P_Bmu]   = CovarianceMatrix[P_Bmu]  [P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMA3];
  CovarianceMatrix[P_Be][P_Bpi]   = CovarianceMatrix[P_Bpi]  [P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMA9];
  CovarianceMatrix[P_Be][P_BK]    = CovarianceMatrix[P_BK]   [P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMA10];
  CovarianceMatrix[P_Be][P_BmuBe] = CovarianceMatrix[P_BmuBe][P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMA3BY5];
  CovarianceMatrix[P_Be][P_BKBpi] = CovarianceMatrix[P_BKBpi][P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMA10BY9];
  CovarianceMatrix[P_Be][P_Bs]    = CovarianceMatrix[P_Bs]   [P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMA110];
  CovarianceMatrix[P_Be][P_Ball]  = CovarianceMatrix[P_Ball] [P_Be] = NodeErrorMatrix[N_GAMMA5][N_GAMMAALL];
  //
  CovarianceMatrix[P_Bmu][P_Bpi]   = CovarianceMatrix[P_Bpi]  [P_Bmu] = NodeErrorMatrix[N_GAMMA3][N_GAMMA9];
  CovarianceMatrix[P_Bmu][P_BK]    = CovarianceMatrix[P_BK]   [P_Bmu] = NodeErrorMatrix[N_GAMMA3][N_GAMMA10];
  CovarianceMatrix[P_Bmu][P_BmuBe] = CovarianceMatrix[P_BmuBe][P_Bmu] = NodeErrorMatrix[N_GAMMA3][N_GAMMA3BY5];
  CovarianceMatrix[P_Bmu][P_BKBpi] = CovarianceMatrix[P_BKBpi][P_Bmu] = NodeErrorMatrix[N_GAMMA3][N_GAMMA10BY9];
  CovarianceMatrix[P_Bmu][P_Bs]    = CovarianceMatrix[P_Bs]   [P_Bmu] = NodeErrorMatrix[N_GAMMA3][N_GAMMA110];
  CovarianceMatrix[P_Bmu][P_Ball]  = CovarianceMatrix[P_Ball] [P_Bmu] = NodeErrorMatrix[N_GAMMA3][N_GAMMAALL];
  //
  CovarianceMatrix[P_Bpi][P_BK]    = CovarianceMatrix[P_BK]   [P_Bpi] = NodeErrorMatrix[N_GAMMA9][N_GAMMA10];
  CovarianceMatrix[P_Bpi][P_BmuBe] = CovarianceMatrix[P_BmuBe][P_Bpi] = NodeErrorMatrix[N_GAMMA9][N_GAMMA3BY5];
  CovarianceMatrix[P_Bpi][P_BKBpi] = CovarianceMatrix[P_BKBpi][P_Bpi] = NodeErrorMatrix[N_GAMMA9][N_GAMMA10BY9];
  CovarianceMatrix[P_Bpi][P_Bs]    = CovarianceMatrix[P_Bs]   [P_Bpi] = NodeErrorMatrix[N_GAMMA9][N_GAMMA110];
  CovarianceMatrix[P_Bpi][P_Ball]  = CovarianceMatrix[P_Ball] [P_Bpi] = NodeErrorMatrix[N_GAMMA9][N_GAMMAALL];
  //
  CovarianceMatrix[P_BK][P_BmuBe]  = CovarianceMatrix[P_BmuBe][P_BK] = NodeErrorMatrix[N_GAMMA10][N_GAMMA3BY5];
  CovarianceMatrix[P_BK][P_BKBpi]  = CovarianceMatrix[P_BKBpi][P_BK] = NodeErrorMatrix[N_GAMMA10][N_GAMMA10BY9];
  CovarianceMatrix[P_BK][P_Bs]     = CovarianceMatrix[P_Bs]   [P_BK] = NodeErrorMatrix[N_GAMMA10][N_GAMMA110];
  CovarianceMatrix[P_BK][P_Ball]   = CovarianceMatrix[P_Ball] [P_BK] = NodeErrorMatrix[N_GAMMA10][N_GAMMAALL];
  //
  CovarianceMatrix[P_BmuBe][P_BKBpi] = CovarianceMatrix[P_BKBpi][P_BmuBe] = NodeErrorMatrix[N_GAMMA3BY5][N_GAMMA10BY9];
  CovarianceMatrix[P_BmuBe][P_Bs]    = CovarianceMatrix[P_Bs]   [P_BmuBe] = NodeErrorMatrix[N_GAMMA3BY5][N_GAMMA110];
  CovarianceMatrix[P_BmuBe][P_Ball]  = CovarianceMatrix[P_Ball] [P_BmuBe] = NodeErrorMatrix[N_GAMMA3BY5][N_GAMMAALL];
  //
  CovarianceMatrix[P_BKBpi][P_Bs]   = CovarianceMatrix[P_Bs]  [P_BKBpi] = NodeErrorMatrix[N_GAMMA10BY9][N_GAMMA110];
  CovarianceMatrix[P_BKBpi][P_Ball] = CovarianceMatrix[P_Ball][P_BKBpi] = NodeErrorMatrix[N_GAMMA10BY9][N_GAMMAALL];
  //
  CovarianceMatrix[P_Bs][P_Ball] = CovarianceMatrix[P_Ball][P_Bs] = NodeErrorMatrix[N_GAMMA110][N_GAMMAALL];
  //
  double** CorrelationMatrix = new double*[N_ParNames]; for (ipar=0;ipar<N_ParNames;++ipar) CorrelationMatrix[ipar] = new double[N_ParNames];
  //
  for (ipar=0;ipar<N_ParNames;++ipar) {
    for (jpar=0;jpar<N_ParNames;++jpar) {
      if (ipar==jpar) {
	CorrelationMatrix[ipar][jpar] = 1;
      } else if (CovarianceMatrix[ipar][ipar] == 0 || CovarianceMatrix[jpar][jpar] == 0) { 
	CorrelationMatrix[ipar][jpar] = 0;
      } else {
	CorrelationMatrix[ipar][jpar] = CovarianceMatrix[ipar][jpar]/sqrt(CovarianceMatrix[ipar][ipar]*CovarianceMatrix[jpar][jpar]);
      }
    }
  }
  //
  // Gamma3by5 : G(mu- nubar(mu) nu(tau)) / G(e- nubar(e) nu(tau))
  //
  fprintf(thisfile,"\n");
  fprintf(thisfile,"%s = %6.4f +- %6.4f\n","B(tau- -> mu- nub nu)/B(tau- -> e- nub nu)", BmuBe, e_BmuBe);
  //
  double fmufe=0;
  double e_fmufe=0;
  double*ecomp_fmufe = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_fmufe, fmufe, e_fmufe, ecomp_fmufe);
  //
  double gmuge=0;
  double e_gmuge=0;
  double*ecomp_gmuge = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_gmuge, gmuge, e_gmuge, ecomp_gmuge);
  fprintf(thisfile,"%s = %8.6f +- %8.6f [Total] +- %8.6f [mtau] +- %8.6f [mmu] +- %8.6f [me]\n",
	  "f(m_mu^2/m_tau^2)/f(m_e^2/m_tau^2)", fmufe, e_fmufe, ecomp_fmufe[P_m_tau], ecomp_fmufe[P_m_mu], ecomp_fmufe[P_m_e]);
  //  for (ipar=0;ipar<N_ParNames;++ipar) {
  //    cout << ipar << " " << CovarianceMatrix[ipar][ipar] << " " << CorrelationMatrix[ipar][ipar]  << " " << ecomp_fmufe[ipar] << " " << endl;
  //  }
  //
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [me]\n",
	  "gmu/ge", gmuge, e_gmuge, ecomp_gmuge[P_BmuBe], ecomp_gmuge[P_m_tau],  ecomp_gmuge[P_m_mu], ecomp_gmuge[P_m_e]);
  //
  // Gamma9: G(tau- --> pi- nu(tau))
  //
  fprintf(thisfile,"\n");
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %%\n","B(tau- -> pi- nu)", 100.*Bpi, 100.*e_Bpi);
  //
  double Bpi_univ = 0;
  double e_Bpi_univ = 0;
  double*ecomp_Bpi_univ = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Bpi_univ, Bpi_univ, e_Bpi_univ, ecomp_Bpi_univ);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [mh] +- %6.3f [tautau] +- %6.3f [tauh] +- %6.3f [Bh] +- %6.3f [Delta]) %%\n",
	  "B(tau- -> pi- nu)_univ", 
	  100.*Bpi_univ, 100.*e_Bpi_univ,
	  100.*ecomp_Bpi_univ[P_m_tau], 
	  100.*ecomp_Bpi_univ[P_m_mu], 
	  100.*ecomp_Bpi_univ[P_m_pim], 
	  100.*ecomp_Bpi_univ[P_tau_tau], 
	  100.*ecomp_Bpi_univ[P_tau_pim], 
	  100.*ecomp_Bpi_univ[P_BR_PimToMumNu], 
	  100.*ecomp_Bpi_univ[P_Delta_TauToPim_over_PimToMu]); 
  //
  double gtaugmu_pi = 0;
  double e_gtaugmu_pi = 0;
  double*ecomp_gtaugmu_pi = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_gtaugmu_pi, gtaugmu_pi, e_gtaugmu_pi, ecomp_gtaugmu_pi);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [mh] +- %6.4f [tautau] +- %6.4f [tauh] +- %6.4f [Bh] +- %6.4f [Delta]\n",
	  "(gtau/gmu)_pi", 
	  gtaugmu_pi, 
	  e_gtaugmu_pi, 
	  ecomp_gtaugmu_pi[P_Bpi],
	  ecomp_gtaugmu_pi[P_m_tau],
	  ecomp_gtaugmu_pi[P_m_mu],
	  ecomp_gtaugmu_pi[P_m_pim],
	  ecomp_gtaugmu_pi[P_tau_tau],
	  ecomp_gtaugmu_pi[P_tau_pim],
	  ecomp_gtaugmu_pi[P_BR_PimToMumNu],
	  ecomp_gtaugmu_pi[P_Delta_TauToPim_over_PimToMu]);
  //
  // Gamma10: G(tau- --> K- nu(tau))
  //
  fprintf(thisfile,"\n");
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %%\n","B(tau- -> K- nu)", 100.*BK, 100.*e_BK);
  //
  double BK_univ = 0;
  double e_BK_univ = 0;
  double*ecomp_BK_univ = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_BK_univ, BK_univ, e_BK_univ, ecomp_BK_univ);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [mh] +- %6.3f [tautau] +- %6.3f [tauh] +- %6.3f [Bh] +- %6.3f [Delta]) %%\n",
	  "B(tau- -> K- nu)_univ", 
	  100.*BK_univ, 100.*e_BK_univ,
	  100.*ecomp_BK_univ[P_m_tau], 
	  100.*ecomp_BK_univ[P_m_mu], 
	  100.*ecomp_BK_univ[P_m_km], 
	  100.*ecomp_BK_univ[P_tau_tau], 
	  100.*ecomp_BK_univ[P_tau_km], 
	  100.*ecomp_BK_univ[P_BR_KmToMumNu], 
	  100.*ecomp_BK_univ[P_Delta_TauToKm_over_KmToMu]); 
  //
  double gtaugmu_K = 0;
  double e_gtaugmu_K = 0;
  double*ecomp_gtaugmu_K = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_gtaugmu_K, gtaugmu_K, e_gtaugmu_K, ecomp_gtaugmu_K);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [mh] +- %6.4f [tautau] +- %6.4f [tauh] +- %6.4f [Bh] +- %6.4f [Delta]\n",
	  "(gtau/gmu)_K",
	  gtaugmu_K,
	  e_gtaugmu_K, 
	  ecomp_gtaugmu_K[P_BK],
	  ecomp_gtaugmu_K[P_m_tau],
	  ecomp_gtaugmu_K[P_m_mu],
	  ecomp_gtaugmu_K[P_m_km],
	  ecomp_gtaugmu_K[P_tau_tau],
	  ecomp_gtaugmu_K[P_tau_km],
	  ecomp_gtaugmu_K[P_BR_KmToMumNu],
	  ecomp_gtaugmu_K[P_Delta_TauToKm_over_KmToMu]);
  //
  // Average of (gtau/gmu)_pi and (gtau/gmu)_K
  //
  fprintf(thisfile,"\n");
  //
  double gtaugmu_pik[2] = {gtaugmu_pi,gtaugmu_K};
  double** cov_gtaugmu_pik = new double*[2]; for (i=0;i<2;++i) cov_gtaugmu_pik[i] = new double[2];
  cov_gtaugmu_pik[0][0] = cov_gtaugmu_pik[1][1] = cov_gtaugmu_pik[0][1] = cov_gtaugmu_pik[1][0] = 0;
  for (ipar=0;ipar<N_ParNames;++ipar) {
    for (jpar=0;jpar<N_ParNames;++jpar) {
      cov_gtaugmu_pik[0][0]+= ecomp_gtaugmu_pi[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_pi[jpar];
      cov_gtaugmu_pik[0][1]+= ecomp_gtaugmu_pi[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_K [jpar];
      cov_gtaugmu_pik[1][0]+= ecomp_gtaugmu_K [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_pi[jpar];
      cov_gtaugmu_pik[1][1]+= ecomp_gtaugmu_K [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_K [jpar];
    }
  }
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","(gtau/gmu)_pi","(gtau/gmu)_K",100.*cov_gtaugmu_pik[0][1]/TMath::Sqrt(cov_gtaugmu_pik[0][0]*cov_gtaugmu_pik[1][1]));
  double gtaugmu_pik_ave = 0;
  double gtaugmu_pik_err = 0;
  double gtaugmu_pik_wt[2] = {0,0};
  calc_average(2,gtaugmu_pik,cov_gtaugmu_pik,gtaugmu_pik_ave,gtaugmu_pik_err,gtaugmu_pik_wt);
  fprintf(thisfile,"%s = %6.4f +- %6.4f has weights = %6.4f [pi], %6.4f [K]\n",
	       "<(gtau/gmu)_pik>",gtaugmu_pik_ave,gtaugmu_pik_err,gtaugmu_pik_wt[0],gtaugmu_pik_wt[1]);
  
  //
  // gtau/gmu = sqrt(Be/Be_from_tautau)
  //
  fprintf(thisfile,"\n");
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %%\n","B(tau- -> e- nub nu)", 100.*Be, 100.*e_Be);
  //
  double Be_from_tautau = 0;
  double e_Be_from_tautau = 0;
  double*ecomp_Be_from_tautau = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Be_from_tautau, Be_from_tautau, e_Be_from_tautau, ecomp_Be_from_tautau);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [me] +- %6.3f [tautau] +- %6.3f [taumu]) %%\n",
	  "B(tau- -> e- nub nu)_tautau ",
	  100.*Be_from_tautau, 
	  100.*e_Be_from_tautau,
	  100.*ecomp_Be_from_tautau[P_m_tau],
	  100.*ecomp_Be_from_tautau[P_m_mu],
	  100.*ecomp_Be_from_tautau[P_m_e],
	  100.*ecomp_Be_from_tautau[P_tau_tau],
	  100.*ecomp_Be_from_tautau[P_tau_mu]);
  //
  double gtaugmu_e = 0;
  double e_gtaugmu_e = 0;
  double*ecomp_gtaugmu_e = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_gtaugmu_e, gtaugmu_e, e_gtaugmu_e, ecomp_gtaugmu_e);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [me] +- %6.4f [tautau] +- %6.4f [taumu] \n",
	  "(gtau/gmu)_e", gtaugmu_e, e_gtaugmu_e, 
	  ecomp_gtaugmu_e[P_Be], 
	  ecomp_gtaugmu_e[P_m_tau], ecomp_gtaugmu_e[P_m_mu], ecomp_gtaugmu_e[P_m_e],
	  ecomp_gtaugmu_e[P_tau_tau], ecomp_gtaugmu_e[P_tau_mu]);
  //
  // Average of (gtau/gmu)_pi and (gtau/gmu)_K and (gtau/gmu)_e
  //
  fprintf(thisfile,"\n");
  //
  double gtaugmu_pike[3] = {gtaugmu_pi,gtaugmu_K,gtaugmu_e};
  double** cov_gtaugmu_pike = new double*[3]; for (i=0;i<3;++i) cov_gtaugmu_pike[i] = new double[3];
  cov_gtaugmu_pike[0][0] = cov_gtaugmu_pike[1][1] = cov_gtaugmu_pike[2][2] = 0;
  cov_gtaugmu_pike[0][1] = cov_gtaugmu_pike[1][0] = 0;
  cov_gtaugmu_pike[0][2] = cov_gtaugmu_pike[2][0] = 0;
  cov_gtaugmu_pike[1][2] = cov_gtaugmu_pike[2][1] = 0;
  for (ipar=0;ipar<N_ParNames;++ipar) {
    for (jpar=0;jpar<N_ParNames;++jpar) {
      cov_gtaugmu_pike[0][0]+= ecomp_gtaugmu_pi[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_pi[jpar];
      cov_gtaugmu_pike[1][1]+= ecomp_gtaugmu_K [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_K [jpar];
      cov_gtaugmu_pike[2][2]+= ecomp_gtaugmu_e [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_e [jpar];
      //
      cov_gtaugmu_pike[0][1]+= ecomp_gtaugmu_pi[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_K [jpar];
      cov_gtaugmu_pike[1][0]+= ecomp_gtaugmu_K [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_pi[jpar];
      //
      cov_gtaugmu_pike[0][2]+= ecomp_gtaugmu_pi[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_e [jpar];
      cov_gtaugmu_pike[2][0]+= ecomp_gtaugmu_e [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_pi[jpar];
      //
      cov_gtaugmu_pike[1][2]+= ecomp_gtaugmu_K [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_e [jpar];
      cov_gtaugmu_pike[2][1]+= ecomp_gtaugmu_e [ipar] * CorrelationMatrix[ipar][jpar] * ecomp_gtaugmu_K [jpar];
    }
  }
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","(gtau/gmu)_pi","(gtau/gmu)_e",100.*cov_gtaugmu_pike[0][2]/TMath::Sqrt(cov_gtaugmu_pike[0][0]*cov_gtaugmu_pike[2][2]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","(gtau/gmu)_K ","(gtau/gmu)_e",100.*cov_gtaugmu_pike[1][2]/TMath::Sqrt(cov_gtaugmu_pike[1][1]*cov_gtaugmu_pike[2][2]));
  //
  double gtaugmu_pike_ave = 0;
  double gtaugmu_pike_err = 0;
  double gtaugmu_pike_wt[3] = {0,0,0};
  calc_average(3,gtaugmu_pike,cov_gtaugmu_pike,gtaugmu_pike_ave,gtaugmu_pike_err,gtaugmu_pike_wt);
  fprintf(thisfile,"%s = %6.4f +- %6.4f has weights = %6.4f [pi], %6.4f [K], %6.4f [e]\n",
	       "<(gtau/gmu)_pike>",gtaugmu_pike_ave,gtaugmu_pike_err,gtaugmu_pike_wt[0],gtaugmu_pike_wt[1],gtaugmu_pike_wt[2]);

  //
  // gtau/ge = sqrt(Bmu/Bmu_from_tautau)
  //
  fprintf(thisfile,"\n");
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %%\n","B(tau- -> mu- nub nu)", 100.*Bmu, 100.*e_Bmu);
  //
  double Bmu_from_tautau = 0;
  double e_Bmu_from_tautau = 0;
  double*ecomp_Bmu_from_tautau = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Bmu_from_tautau, Bmu_from_tautau, e_Bmu_from_tautau, ecomp_Bmu_from_tautau);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [me] +- %6.3f [tautau] +- %6.3f [taumu]) %%\n",
	  "B(tau- -> mu- nub nu)_tautau ",
	  100.*Bmu_from_tautau, 
	  100.*e_Bmu_from_tautau,
	  100.*ecomp_Bmu_from_tautau[P_m_tau],
	  100.*ecomp_Bmu_from_tautau[P_m_mu],
	  100.*ecomp_Bmu_from_tautau[P_m_e],
	  100.*ecomp_Bmu_from_tautau[P_tau_tau],
	  100.*ecomp_Bmu_from_tautau[P_tau_mu]);
  //
  double gtauge_mu = 0;
  double e_gtauge_mu = 0;
  double*ecomp_gtauge_mu = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_gtauge_mu, gtauge_mu, e_gtauge_mu, ecomp_gtauge_mu);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [me] +- %6.4f [tautau] +- %6.4f [taumu] \n",
	  "(gtau/ge)_mu", gtauge_mu, e_gtauge_mu, 
	  ecomp_gtauge_mu[P_Bmu], 
	  ecomp_gtauge_mu[P_m_tau], ecomp_gtauge_mu[P_m_mu], ecomp_gtauge_mu[P_m_e],
	  ecomp_gtauge_mu[P_tau_tau], ecomp_gtauge_mu[P_tau_mu]);
  //
  // Gamma10: G(tau- --> K- nu(tau))
  //
  fprintf(thisfile,"\n");
  //
  double   Vus_TauToKmNu = 0;
  double e_Vus_TauToKmNu = 0;
  double *ecomp_Vus_TauToKmNu = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Vus_K, Vus_TauToKmNu, e_Vus_TauToKmNu, ecomp_Vus_TauToKmNu);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [tautau] +- %6.4f [mK] +- %6.4f [fK]\n",
	  "|Vus|_TauToKmNu",
	  Vus_TauToKmNu,
	  e_Vus_TauToKmNu,
	  ecomp_Vus_TauToKmNu[P_BK],
	  ecomp_Vus_TauToKmNu[P_m_tau],
	  ecomp_Vus_TauToKmNu[P_tau_tau],
	  ecomp_Vus_TauToKmNu[P_m_km],
	  ecomp_Vus_TauToKmNu[P_fK]);
  //
  // Gamma10by9 : G(K- nu(tau)) / G(pi- nu(tau))
  //
  fprintf(thisfile,"\n");
  fprintf(thisfile, "%s = %6.4f +- %6.4f\n","G(K- nu(tau)) / G(pi- nu(tau))", BKBpi, e_BKBpi);
  //
  double   Vus_TauToKmOverPimNu = 0;
  double e_Vus_TauToKmOverPimNu = 0;
  double *ecomp_Vus_TauToKmOverPimNu = new double[N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Vus_Kpi, Vus_TauToKmOverPimNu, e_Vus_TauToKmOverPimNu, ecomp_Vus_TauToKmOverPimNu);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mK] +- %6.4f [mpi] +- %6.4f [fK/fpi] +- %6.4f [Vud] +- %6.4f [Delta]\n",
	  "|Vus|_TauToK/Pi",
	  Vus_TauToKmOverPimNu,
	  e_Vus_TauToKmOverPimNu,
	  ecomp_Vus_TauToKmOverPimNu[P_BKBpi],
	  ecomp_Vus_TauToKmOverPimNu[P_m_tau],
          ecomp_Vus_TauToKmOverPimNu[P_m_km],
          ecomp_Vus_TauToKmOverPimNu[P_m_pim],
	  ecomp_Vus_TauToKmOverPimNu[P_fKfpi],
	  ecomp_Vus_TauToKmOverPimNu[P_Vud],
	  ecomp_Vus_TauToKmOverPimNu[P_Delta_Kpi]);
  //
  // Be from Bmu
  //
  fprintf(thisfile,"\n");
  //
  double Be_from_Bmu = 0;
  double e_Be_from_Bmu = 0;
  double*ecomp_Be_from_Bmu = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Be_from_Bmu, Be_from_Bmu, e_Be_from_Bmu, ecomp_Be_from_Bmu);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [Btau] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [me]) %%\n",
	  "B(tau- -> e- nub nu)_Bmu ",
	  100.*Be_from_Bmu, 100.*e_Be_from_Bmu, 100.*ecomp_Be_from_Bmu[P_Bmu], 100.*ecomp_Be_from_Bmu[P_m_tau], 100.*ecomp_Be_from_Bmu[P_m_mu], 100.*ecomp_Be_from_Bmu[P_m_e]);
  //
  // Be from Bpi
  //
  fprintf(thisfile,"\n");
  double Be_from_Bpi=0;
  double e_Be_from_Bpi=0;
  double *ecomp_Be_from_Bpi = new double[N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Be_from_Bpi, Be_from_Bpi, e_Be_from_Bpi, ecomp_Be_from_Bpi);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mh] +- %6.3f [tauh] +- %6.3f [Btau] +- %6.3f [Bh] +- %6.3f [Delta]) %%\n",
	  "B(tau- -> e- nub nu)_Bpi ",
	  100.*Be_from_Bpi, 100.*e_Be_from_Bpi, 
	  100.*ecomp_Be_from_Bpi[P_m_tau], 100.*ecomp_Be_from_Bpi[P_m_pim], 100.*ecomp_Be_from_Bpi[P_tau_pim], 
	  100.*ecomp_Be_from_Bpi[P_Bpi], 100.*ecomp_Be_from_Bpi[P_BR_PimToMumNu], 100.*ecomp_Be_from_Bpi[P_Delta_TauToPim_over_PimToMu]);
  //
  // Be from BK
  //
  fprintf(thisfile,"\n");
  //
  double Be_from_BK=0;
  double e_Be_from_BK=0;
  double *ecomp_Be_from_BK = new double[N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Be_from_BK, Be_from_BK, e_Be_from_BK, ecomp_Be_from_BK);
  fprintf(thisfile,"%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mh] +- %6.3f [tauh] +- %6.3f [Btau] +- %6.3f [Bh] +- %6.3f [Delta]) %%\n",
	  "B(tau- -> e- nub nu)_BK  ",
	  100.*Be_from_BK, 100.*e_Be_from_BK, 
	  100.*ecomp_Be_from_BK[P_m_tau], 100.*ecomp_Be_from_BK[P_m_km], 100.*ecomp_Be_from_BK[P_tau_km], 
	  100.*ecomp_Be_from_BK[P_BK], 100.*ecomp_Be_from_BK[P_BR_KmToMumNu], 100.*ecomp_Be_from_BK[P_Delta_TauToKm_over_KmToMu]);

  //
  // Average of Be and Be_from_Bmu and Be_from_tautau and Bpi and BK
  //
  fprintf(thisfile,"\n");
  double*ecomp_Be = new double [N_ParNames];
  for (ipar=0;ipar<N_ParNames;++ipar) {ecomp_Be[ipar] = 0;} ; ecomp_Be[P_Be] = e_Be;
  //
  double Be_univ_val[5] = {Be, Be_from_Bmu, Be_from_tautau, Be_from_Bpi,  Be_from_BK};
  double** Be_univ_cov = new double*[5]; for (i=0;i<5;++i) Be_univ_cov[i] = new double[5];
  Be_univ_cov[0][0] = Be_univ_cov[1][1] = Be_univ_cov[2][2] = Be_univ_cov[3][3] = Be_univ_cov[4][4] = 0;
  Be_univ_cov[0][1] = Be_univ_cov[1][0] = 0;
  Be_univ_cov[0][2] = Be_univ_cov[2][0] = 0;
  Be_univ_cov[0][3] = Be_univ_cov[3][0] = 0;
  Be_univ_cov[0][4] = Be_univ_cov[4][0] = 0;
  Be_univ_cov[1][2] = Be_univ_cov[2][1] = 0;
  Be_univ_cov[1][3] = Be_univ_cov[3][1] = 0;
  Be_univ_cov[1][4] = Be_univ_cov[4][1] = 0;
  Be_univ_cov[2][3] = Be_univ_cov[3][2] = 0;
  Be_univ_cov[2][4] = Be_univ_cov[4][2] = 0;
  Be_univ_cov[3][4] = Be_univ_cov[4][3] = 0;
  //
  for (ipar=0;ipar<N_ParNames;++ipar) {
    for (jpar=0;jpar<N_ParNames;++jpar) {
      Be_univ_cov[0][0] += ecomp_Be[ipar]             * CorrelationMatrix[ipar][jpar] * ecomp_Be[jpar]; 
      Be_univ_cov[1][1] += ecomp_Be_from_Bmu[ipar]    * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_Bmu[jpar] ;
      Be_univ_cov[2][2] += ecomp_Be_from_tautau[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_tautau[jpar] ;
      Be_univ_cov[3][3] += ecomp_Be_from_Bpi[ipar]    * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_Bpi[jpar] ;
      Be_univ_cov[4][4] += ecomp_Be_from_BK[ipar]     * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_BK[jpar] ;
      //
      Be_univ_cov[0][1] += ecomp_Be[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_Bmu[jpar] ;
      Be_univ_cov[0][2] += ecomp_Be[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_tautau[jpar] ;
      Be_univ_cov[0][3] += ecomp_Be[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_Bpi[jpar] ;
      Be_univ_cov[0][4] += ecomp_Be[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_BK[jpar] ;
      //
      Be_univ_cov[1][2] += ecomp_Be_from_Bmu[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_tautau[jpar] ;
      Be_univ_cov[1][3] += ecomp_Be_from_Bmu[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_Bpi[jpar] ;
      Be_univ_cov[1][4] += ecomp_Be_from_Bmu[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_BK[jpar] ;
      //
      Be_univ_cov[2][3] += ecomp_Be_from_tautau[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_Bpi[jpar] ;
      Be_univ_cov[2][4] += ecomp_Be_from_tautau[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_BK[jpar] ;
      //
      Be_univ_cov[3][4] += ecomp_Be_from_Bpi[ipar] * CorrelationMatrix[ipar][jpar] * ecomp_Be_from_BK[jpar] ;
    }
    /*
    cout << " ipar = " << ipar 
	 << " Be_univ_cov[1][3] = " << 100.*Be_univ_cov[1][3] / (e_Be_from_Bpi * e_Be_from_Bmu) << " % "
	 << " Be_univ_cov[3][4] = " << 100.*Be_univ_cov[3][4] / (e_Be_from_Bpi * e_Be_from_BK)  << " % "
	 << endl;
    */
  }
  Be_univ_cov[1][0] = Be_univ_cov[0][1] ;
  Be_univ_cov[2][0] = Be_univ_cov[0][2] ;
  Be_univ_cov[3][0] = Be_univ_cov[0][3] ;
  Be_univ_cov[4][0] = Be_univ_cov[0][4] ;
  //
  Be_univ_cov[2][1] = Be_univ_cov[1][2] ;
  Be_univ_cov[3][1] = Be_univ_cov[1][3] ;
  Be_univ_cov[4][1] = Be_univ_cov[1][4] ;
  //
  Be_univ_cov[3][2] = Be_univ_cov[2][3] ;
  Be_univ_cov[4][2] = Be_univ_cov[2][4] ;
  //
  Be_univ_cov[4][3] = Be_univ_cov[3][4] ;
  //
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_Bmu","Be            ",100.*Be_univ_cov[0][1]/TMath::Sqrt(Be_univ_cov[0][0]*Be_univ_cov[1][1]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_Bmu","Be_from_tautau",100.*Be_univ_cov[2][1]/TMath::Sqrt(Be_univ_cov[2][2]*Be_univ_cov[1][1]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_Bpi","Be            ",100.*Be_univ_cov[0][3]/TMath::Sqrt(Be_univ_cov[0][0]*Be_univ_cov[3][3]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_Bpi","Be_from_Bmu   ",100.*Be_univ_cov[1][3]/TMath::Sqrt(Be_univ_cov[1][1]*Be_univ_cov[3][3]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_Bpi","Be_from_tautau",100.*Be_univ_cov[2][3]/TMath::Sqrt(Be_univ_cov[2][2]*Be_univ_cov[3][3]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_BK ","Be            ",100.*Be_univ_cov[0][4]/TMath::Sqrt(Be_univ_cov[0][0]*Be_univ_cov[4][4]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_BK ","Be_from_Bmu   ",100.*Be_univ_cov[1][4]/TMath::Sqrt(Be_univ_cov[1][1]*Be_univ_cov[4][4]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_BK ","Be_from_tautau",100.*Be_univ_cov[2][4]/TMath::Sqrt(Be_univ_cov[2][2]*Be_univ_cov[4][4]));
  fprintf(thisfile,"Corr between %s and %s = %6.2f %%\n","Be_from_BK ","Be_from_Bpi   ",100.*Be_univ_cov[3][4]/TMath::Sqrt(Be_univ_cov[3][3]*Be_univ_cov[4][4]));
  //
  double Be_univ3 = 0;
  double Be_univ3_err = 0;
  double Be_univ3_wt[3] = {0,0,0};
  calc_average(3,Be_univ_val,Be_univ_cov,Be_univ3,Be_univ3_err,Be_univ3_wt);
  //
  // Cross-check the error on Be_univ3 assuming the weights have no errors
  //    f = a A + b B + c C => e_f^2 = a^2 e_A^2 + b^2 e_B^2 + c^2 e_C^2 
  //                                 + 2 ab Corr(A,B) e_A e_B + 2 bc Corr(B,C) e_B e_C + 2 ca Corr(C,A) e_C e_A
  // N.B.:  The weights are function of errors, eg. e[mtau], e[tautau], e[Be], e[Bmu], etc.; so it is reasonable assumption that error on these errors = 0.
  //                                             
  double Be_univ3_err_re = TMath::Sqrt( TMath::Power(Be_univ3_wt[0],2) * Be_univ_cov[0][0] + 
					TMath::Power(Be_univ3_wt[1],2) * Be_univ_cov[1][1] + 
					TMath::Power(Be_univ3_wt[2],2) * Be_univ_cov[2][2] + 
					2 * Be_univ3_wt[0] * Be_univ3_wt[1] * Be_univ_cov[0][1] + 
					2 * Be_univ3_wt[0] * Be_univ3_wt[2] * Be_univ_cov[0][2] + 
					2 * Be_univ3_wt[2] * Be_univ3_wt[1] * Be_univ_cov[2][1] );
  if (fabs(Be_univ3_err - Be_univ3_err_re)>.001) {
    cout << " Be_univ3_err =  " << Be_univ3_err << " Be_univ3_err_re = " << Be_univ3_err_re << endl;
    abort();
  }
  //
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %% has weights = %6.4f [Be], %6.4f [Bmu], %6.4f [tautau]\n",
	  "<B(tau- -> e- nub nu)_univ3>",100.*Be_univ3,100.*Be_univ3_err, Be_univ3_wt[0], Be_univ3_wt[1], Be_univ3_wt[2]);
  //
  double Be_univ4 = 0;
  double Be_univ4_err = 0;
  double Be_univ4_wt[4] = {0,0,0,0};
  calc_average(4,Be_univ_val,Be_univ_cov,Be_univ4,Be_univ4_err,Be_univ4_wt);
  //
  // Cross-check the error on Be_univ4 assuming the weights have no errors
  //    f = a A + b B + c C + d D => e_f^2 = a^2 e_A^2 + b^2 e_B^2 + c^2 e_C^2 + d^2 e_D^2 
  //                                       + 2 ab Corr(A,B) e_A e_B + 2 bc Corr(B,C) e_B e_C + 2 ca Corr(C,A) e_C e_A
  //                                       + 2 ad Corr(A,D) e_A e_D + 2 bd Corr(B,D) e_B e_D + 2 cd Corr(C,D) e_C e_D
  // N.B.:  The weights are function of errors, eg. e[mtau], e[tautau], e[Be], e[Bmu], etc.; so it is reasonable assumption that error on these errors = 0.
  //                                             
  double Be_univ4_err_re = TMath::Sqrt( TMath::Power(Be_univ4_wt[0],2) * Be_univ_cov[0][0] + 
					TMath::Power(Be_univ4_wt[1],2) * Be_univ_cov[1][1] + 
					TMath::Power(Be_univ4_wt[2],2) * Be_univ_cov[2][2] + 
					TMath::Power(Be_univ4_wt[3],2) * Be_univ_cov[3][3] + 
					2 * Be_univ4_wt[0] * Be_univ4_wt[1] * Be_univ_cov[0][1] + 
					2 * Be_univ4_wt[0] * Be_univ4_wt[2] * Be_univ_cov[0][2] + 
					2 * Be_univ4_wt[0] * Be_univ4_wt[3] * Be_univ_cov[0][3] + 
					2 * Be_univ4_wt[1] * Be_univ4_wt[2] * Be_univ_cov[1][2] +
					2 * Be_univ4_wt[1] * Be_univ4_wt[3] * Be_univ_cov[1][3] + 
					2 * Be_univ4_wt[2] * Be_univ4_wt[3] * Be_univ_cov[2][3] );
  if (fabs(Be_univ4_err - Be_univ4_err_re)>0.001) {
    cout << " Be_univ4_err = " << Be_univ4_err  << " Be_univ4_err_re = " << Be_univ4_err_re << endl;
    abort();
  }
  //
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %% has weights = %6.4f [Be], %6.4f [Bmu], %6.4f [tautau], %6.4f [Bpi]\n",  
	  "<B(tau- -> e- nub nu)_univ4>",100.*Be_univ4,100.*Be_univ4_err, Be_univ4_wt[0], Be_univ4_wt[1], Be_univ4_wt[2], Be_univ4_wt[3]);
  //
  double Be_univ5 = 0;
  double Be_univ5_err = 0;
  double Be_univ5_wt[5] = {0,0,0,0,0};
  calc_average(5,Be_univ_val,Be_univ_cov,Be_univ5,Be_univ5_err,Be_univ5_wt);
  //
  // Cross-check the error on Be_univ5 assuming the weights have no errors
  //    f = a A + b B + c C + d D + e E => e_f^2 = a^2 e_A^2 + b^2 e_B^2 + c^2 e_C^2 + d^2 e_D^2 + e^2 e_E^2 
  //                                       + 2 ab Corr(A,B) e_A e_B + 2 bc Corr(B,C) e_B e_C + 2 ca Corr(C,A) e_C e_A
  //                                       + 2 ad Corr(A,D) e_A e_D + 2 bd Corr(B,D) e_B e_D + 2 cd Corr(C,D) e_C e_D
  //                                       + 2 ae Corr(A,E) e_A e_E + 2 be Corr(B,E) e_B e_E + 2 ce Corr(C,E) e_C e_E + 2 de Corr(D,E) e_D e_E
  // N.B.:  The weights are function of errors, eg. e[mtau], e[tautau], e[Be], e[Bmu], etc.; so it is reasonable assumption that error on these errors = 0.
  //                                             
  double Be_univ5_err_re = TMath::Sqrt( TMath::Power(Be_univ5_wt[0],2) * Be_univ_cov[0][0] + 
					TMath::Power(Be_univ5_wt[1],2) * Be_univ_cov[1][1] + 
					TMath::Power(Be_univ5_wt[2],2) * Be_univ_cov[2][2] + 
					TMath::Power(Be_univ5_wt[3],2) * Be_univ_cov[3][3] + 
					TMath::Power(Be_univ5_wt[4],2) * Be_univ_cov[4][4] + 
					2 * Be_univ5_wt[0] * Be_univ5_wt[1] * Be_univ_cov[0][1] + 
					2 * Be_univ5_wt[0] * Be_univ5_wt[2] * Be_univ_cov[0][2] + 
					2 * Be_univ5_wt[0] * Be_univ5_wt[3] * Be_univ_cov[0][3] + 
					2 * Be_univ5_wt[0] * Be_univ5_wt[4] * Be_univ_cov[0][4] +
					2 * Be_univ5_wt[1] * Be_univ5_wt[2] * Be_univ_cov[1][2] +
					2 * Be_univ5_wt[1] * Be_univ5_wt[3] * Be_univ_cov[1][3] + 
					2 * Be_univ5_wt[1] * Be_univ5_wt[4] * Be_univ_cov[1][4] +
					2 * Be_univ5_wt[2] * Be_univ5_wt[3] * Be_univ_cov[2][3] +
					2 * Be_univ5_wt[2] * Be_univ5_wt[4] * Be_univ_cov[2][4] +
					2 * Be_univ5_wt[3] * Be_univ5_wt[4] * Be_univ_cov[3][4] );
  if (fabs(Be_univ5_err - Be_univ5_err_re)>0.001) {
    cout << "Be_univ5_err = " << Be_univ5_err << " Be_univ5_err_re = " << Be_univ5_err_re << endl;
    abort();
  }
  //
  fprintf(thisfile,"%s = (%6.3f +- %6.3f) %% has weights = %6.4f [Be], %6.4f [Bmu], %6.4f [tautau], %6.4f [Bpi], %6.4f [BK]\n", 
	  "<B(tau- -> e- nub nu)_univ5>", 100.*Be_univ5,100.*Be_univ5_err, Be_univ5_wt[0], Be_univ5_wt[1], Be_univ5_wt[2], Be_univ5_wt[3], Be_univ5_wt[4]) ;
  //
  par[P_wt_Be] = Be_univ5_wt[0];
  par[P_wt_Be_from_Bmu] = Be_univ5_wt[1];
  par[P_wt_Be_from_tautau] = Be_univ5_wt[2];
  par[P_wt_Be_from_Bpi] = Be_univ5_wt[3];
  par[P_wt_Be_from_BK] = Be_univ5_wt[4];
  //
  // Bhadrons
  //
  fprintf(thisfile,"\n");
  double   Bhadrons = 0;
  double e_Bhadrons = 0;
  double*ecomp_Bhadrons = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Bhadrons, Bhadrons, e_Bhadrons, ecomp_Bhadrons);
  fprintf(thisfile,
	  "%s = (%6.3f +- %6.3f [Total] +- %6.3f [tautau] +- %6.3f [mtau] +- %6.3f [Be] +- %6.3f [Bmu] +- %6.3f [Bpi] +- %6.3f [BK] +- %6.3f [Bs] +- %6.3f [Ball]) %%\n",
  	  "B(tau -> hadrons)", 100.*Bhadrons, 100.*e_Bhadrons,
  	  100.*ecomp_Bhadrons[P_tau_tau], 100.*ecomp_Bhadrons[P_m_tau], 
	  100.*ecomp_Bhadrons[P_Be], 100.*ecomp_Bhadrons[P_Bmu], 100.*ecomp_Bhadrons[P_Bpi], 100.*ecomp_Bhadrons[P_BK], 100.*ecomp_Bhadrons[P_Bs], 100.*ecomp_Bhadrons[P_Ball]);
  //
  // Rhadrons
  //
  double   Rhadrons = 0;
  double e_Rhadrons = 0;
  double*ecomp_Rhadrons = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Rhadrons, Rhadrons, e_Rhadrons, ecomp_Rhadrons);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [tautau] +- %6.4f [mtau] +- %6.4f [Be] +- %6.4f [Bmu] +- %6.4f [Bpi] +- %6.4f [BK] +- %6.4f [Bs] +- %6.4f [Ball]\n",
	  "R(tau -> hadrons)", Rhadrons, e_Rhadrons,
  	  ecomp_Rhadrons[P_tau_tau], ecomp_Rhadrons[P_m_tau], 
	  ecomp_Rhadrons[P_Be], ecomp_Rhadrons[P_Bmu], ecomp_Rhadrons[P_Bpi], ecomp_Rhadrons[P_BK], ecomp_Rhadrons[P_Bs], ecomp_Rhadrons[P_Ball]);
  //
  // Rstrange
  //

  fprintf(thisfile,"\n%s = (%6.3f +- %6.3f) %% \n","B(tau -> strange)", 100.*Bstrange, 100.*e_Bstrange);
  double Rstrange = 0;
  double e_Rstrange = 0;
  double*ecomp_Rstrange = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Rstrange, Rstrange, e_Rstrange, ecomp_Rstrange);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [tautau] +- %6.4f [mtau] +- %6.4f [Be] +- %6.4f [Bmu] +- %6.4f [Bpi] +- %6.4f [BK] +- %6.4f [Bs] + %6.4f [Ball]\n",
	  "R(tau -> strange)", Rstrange, e_Rstrange,
	  ecomp_Rstrange[P_tau_tau], ecomp_Rstrange[P_m_tau], 
	  ecomp_Rstrange[P_Be], ecomp_Rstrange[P_Bmu], ecomp_Rstrange[P_Bpi], ecomp_Rstrange[P_BK], ecomp_Rstrange[P_Bs], ecomp_Rstrange[P_Ball]);
  //
  // Rnonstrange
  //
  double Rnonstrange = 0;
  double e_Rnonstrange = 0;
  double*ecomp_Rnonstrange = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Rnonstrange, Rnonstrange, e_Rnonstrange, ecomp_Rnonstrange);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [tautau] +- %6.4f [mtau] +- %6.4f [Be] +- %6.4f [Bmu] +- %6.4f [Bpi] +- %6.4f [BK] +- %6.4f [Bs] +- %6.4f [Ball]\n",
	  "R(tau -> non-str)", Rnonstrange, e_Rnonstrange,
	  ecomp_Rnonstrange[P_tau_tau], ecomp_Rnonstrange[P_m_tau], 
	  ecomp_Rnonstrange[P_Be], ecomp_Rnonstrange[P_Bmu], ecomp_Rnonstrange[P_Bpi], ecomp_Rnonstrange[P_BK], ecomp_Rnonstrange[P_Bs], ecomp_Rnonstrange[P_Ball]);
  //
  // Vus_strange
  //
  fprintf(thisfile,"\n%s = (%6.3f +- %6.3f) %% \n","B(tau -> total)", 100.*Btotal, 100.*e_Btotal);
  double Vus_strange = 0;
  double e_Vus_strange = 0;
  double*ecomp_Vus_strange = new double [N_ParNames];
  calc_generic(par, epar, CovarianceMatrix, func_Vus_strange, Vus_strange, e_Vus_strange, ecomp_Vus_strange);
  fprintf(thisfile,"%s = %6.4f +- %6.4f [Total] +- %6.4f [tautau] +- %6.4f [mtau] +- %6.4f [Be] +- %6.4f [Bmu] +- %6.4f [Bpi] +- %6.4f [BK] +- %6.4f [Bs] +- %6.4f [Ball] +- %6.4f [Vud] +- %6.4f [ms]\n",
	  "|Vus|_strange", Vus_strange, e_Vus_strange,
	  ecomp_Vus_strange[P_tau_tau], ecomp_Vus_strange[P_m_tau], 
	  ecomp_Vus_strange[P_Be], ecomp_Vus_strange[P_Bmu], ecomp_Vus_strange[P_Bpi], ecomp_Vus_strange[P_BK], ecomp_Vus_strange[P_Bs], ecomp_Vus_strange[P_Ball],
	  ecomp_Vus_strange[P_Vud], ecomp_Vus_strange[P_Delta_Rth]);
  fprintf(thisfile,"Relative Error (in %%): +- %5.2f [Total] +- %5.2f [tautau] +- %5.2f [mtau] +- %5.2f [Be] +- %5.2f [Bmu] +- %6.4f [Bpi] +- %6.4f [BK] +- %5.2f [Bs] +- %5.2f [Ball] +- %5.2f [Vud] +- %5.2f [ms]\n", 
	  e_Vus_strange*100./Vus_strange,
	  ecomp_Vus_strange[P_tau_tau]*100./Vus_strange, ecomp_Vus_strange[P_m_tau]*100./Vus_strange, 
	  ecomp_Vus_strange[P_Be]*100./Vus_strange, ecomp_Vus_strange[P_Bmu]*100./Vus_strange, 
	  ecomp_Vus_strange[P_Bpi]*100./Vus_strange, ecomp_Vus_strange[P_BK]*100./Vus_strange, 
	  ecomp_Vus_strange[P_Bs]*100./Vus_strange, ecomp_Vus_strange[P_Ball]*100./Vus_strange, 
	  ecomp_Vus_strange[P_Vud]*100./Vus_strange, ecomp_Vus_strange[P_Delta_Rth]*100./Vus_strange);
  //
  // Vus_unitarity
  //
  fprintf(thisfile,"\n");
  double Vus_unitarity = TMath::Sqrt(1-Vud*Vud);
  double e_Vus_unitarity = -e_Vud*Vud/TMath::Sqrt(1-Vud*Vud);
  //
  double diff_Vus_TauToKmNu = Vus_TauToKmNu - Vus_unitarity;
  double e_diff_Vus_TauToKmNu = TMath::Sqrt(e_Vus_TauToKmNu * e_Vus_TauToKmNu + e_Vus_unitarity * e_Vus_unitarity);
  //
  double diff_Vus_TauToKmOverPimNu = Vus_TauToKmOverPimNu - Vus_unitarity;
  double e_diff_Vus_TauToKmOverPimNu = TMath::Sqrt(e_Vus_TauToKmOverPimNu * e_Vus_TauToKmOverPimNu + e_Vus_unitarity * e_Vus_unitarity - 2 * ecomp_Vus_TauToKmOverPimNu[P_Vud] * e_Vus_unitarity);
  //
  double diff_Vus_strange = Vus_strange - Vus_unitarity;
  double e_diff_Vus_strange = TMath::Sqrt(e_Vus_strange * e_Vus_strange + e_Vus_unitarity * e_Vus_unitarity - 2 * ecomp_Vus_strange[P_Vud] * e_Vus_unitarity);
  //
  fprintf(thisfile,"|Vus_unitarity| = %6.4f +- %6.4f [Vud]; Difference w.r.t unitarity: |Vus|_TauToKmNu = %3.1f |Vus|_TauToK/Pi = %3.1f |Vus|_strange = %3.1f \n",
	  Vus_unitarity, e_Vus_unitarity, 
	  diff_Vus_TauToKmNu/e_diff_Vus_TauToKmNu, 
	  diff_Vus_TauToKmOverPimNu/e_diff_Vus_TauToKmOverPimNu,
	  diff_Vus_strange/e_diff_Vus_strange);
  //
  // Clean Up
  //
  delete [] ecomp_fmufe;
  delete [] ecomp_gmuge;
  delete [] ecomp_Be_from_Bmu;
  delete [] ecomp_Be_from_Bpi;
  delete [] ecomp_Be_from_BK;
  delete [] ecomp_Bpi_univ;
  delete [] ecomp_BK_univ;
  delete [] ecomp_gtaugmu_pi;
  delete [] ecomp_gtaugmu_K;
  delete [] ecomp_Be_from_tautau;
  delete [] ecomp_Bmu_from_tautau;
  delete [] ecomp_gtauge_mu;
  delete [] ecomp_Vus_TauToKmNu; 
  delete [] ecomp_Vus_TauToKmOverPimNu;
  delete [] ecomp_Be;
  delete [] ecomp_Bhadrons;
  delete [] ecomp_Rhadrons;
  delete [] ecomp_Rstrange;
  delete [] ecomp_Rnonstrange;
  delete [] ecomp_Vus_strange;
  delete [] CorrelationMatrix;
  delete [] CovarianceMatrix;
  delete [] epar;
  delete [] par;
  delete [] Be_univ_cov;
  delete [] cov_gtaugmu_pike;
  delete [] cov_gtaugmu_pik;
}
// ----------------------------------------------------------------------
void print_node_def(int nnode, char** a_nodename, string* nodetitle,  vector<string> nodegammaname, bool* node_is_base, 
		    int ncoef, double* coef, string* coefnames,
		    int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		    int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef, 
		    vector<int> baseparm, vector<int> basegamma, string* basetitle){
  //
  int p, inode, ipar, icoef;
  FILE *nodefile[3];
  nodefile[0]=fopen("all_node_def.txt","w");
  nodefile[1]=fopen("derived_node_def.txt","w");
  nodefile[2]=fopen("base_node_def.txt","w");
  for (p=0;p<3;++p) {
    for (inode=0;inode<nnode;++inode) {
      if (node_is_base[inode]) {
	if (p==1) continue;
      } else {
	if (p==2) continue;
      }
      // Print a commented line about this node
      fprintf (nodefile[p], "\n* %-7s : %s : %s\n",a_nodename[inode],nodetitle[inode].data(), nodegammaname[inode].data());
      //
      // Print node expressed in terms of base decay mode names
      //
      fprintf (nodefile[p], "\n  %s =\n", nodetitle[inode].data());
      for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
        if (ipar==0) fprintf (nodefile[p], "("); else fprintf (nodefile[p], " ");
        int parm=node_num_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
	for (icoef=0;icoef<ncoef;++icoef) if (coef[icoef] == node_num_coef[inode].at(ipar) ) break;
	fprintf (nodefile[p], " G(%s) * %s ", basetitle[quan-1].data(), coefnames[icoef].data());
        if (ipar==node_num_parm[inode].size()-1) fprintf (nodefile[p], " )"); else fprintf (nodefile[p], " +\n");
      }
      if (node_den_parm[inode].size()==0) fprintf (nodefile[p], "\n"); else fprintf (nodefile[p], " /\n"); 
      for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
        if (ipar==0) fprintf (nodefile[p], "("); else fprintf (nodefile[p], " ");
        int parm=node_den_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
	for (icoef=0;icoef<ncoef;++icoef) if (coef[icoef] == node_den_coef[inode].at(ipar) ) break;
        fprintf (nodefile[p], " G(%s) * %s", basetitle[quan-1].data(), coefnames[icoef].data());
        if (ipar==node_den_parm[inode].size()-1) fprintf (nodefile[p], " )\n"); else fprintf (nodefile[p], " +\n");
      }
      //
      // Print node expressed in terms of base gamma numbers
      //
      fprintf (nodefile[p], "\n%s = ", nodegammaname[inode].data());
      for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
        if (ipar==0) { fprintf (nodefile[p], "(") ; } else {fprintf (nodefile[p], " + ") ;}
        int parm=node_num_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
        fprintf (nodefile[p], "Gamma%d*%12.10f", basegamma[quan-1], node_num_coef[inode].at(ipar));
        if (ipar==node_num_parm[inode].size()-1) fprintf (nodefile[p], ")");
      }
      if (node_den_parm[inode].size()==0) fprintf (nodefile[p], "\n") ; 
      for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
        if (ipar==0) { fprintf (nodefile[p], " / (") ; } else {fprintf (nodefile[p], " + ") ;}
        int parm=node_den_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
        fprintf (nodefile[p], "Gamma%d*%12.10f", basegamma[quan-1], node_den_coef[inode].at(ipar));
        if (ipar==node_den_parm[inode].size()-1) fprintf (nodefile[p], ")\n");
      }
      //
      // Print node expressed in terms of base parameter numbers
      //
      fprintf (nodefile[p], "\n%s = ", a_nodename[inode]);
      for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
        if (ipar==0) { fprintf (nodefile[p], "(") ; } else {fprintf (nodefile[p], " + ") ;}
        int parm=node_num_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
        fprintf (nodefile[p], "Parm%d*%12.10f", baseparm[quan-1], node_num_coef[inode].at(ipar));
        if (ipar==node_num_parm[inode].size()-1) fprintf (nodefile[p], ")");
      }
      if (node_den_parm[inode].size()==0) fprintf (nodefile[p], "\n") ; 
      for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
        if (ipar==0) { fprintf (nodefile[p], " / (") ; } else {fprintf (nodefile[p], " + ") ;}
        int parm=node_den_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
        fprintf (nodefile[p], "Parm%d*%12.10f", baseparm[quan-1], node_den_coef[inode].at(ipar));
        if (ipar==node_den_parm[inode].size()-1) fprintf (nodefile[p], ")\n");
      }
      //
    }
    fclose(nodefile[p]);
  }
}
// ----------------------------------------------------------------------
void define_nodes(int nnode, vector<int> * node_num_parm, vector<int> * node_den_parm, vector<int> baseparm,
		  vector<int> * node_parm, // output
		  vector<int> * node_quan, // output
		  int * first_quan // output
		  ){
  int inode, ipar;
  for (inode=0;inode<nnode;++inode){
    node_parm[inode].insert(node_parm[inode].end(),node_num_parm[inode].begin(),node_num_parm[inode].end());
    node_parm[inode].insert(node_parm[inode].end(),node_den_parm[inode].begin(),node_den_parm[inode].end());
    sort(node_parm[inode].begin(), node_parm[inode].end());//sort must be called before unique
    vector<int>::iterator new_end=unique(node_parm[inode].begin(), node_parm[inode].end());
    node_parm[inode].erase(new_end, node_parm[inode].end()); // <--
    first_quan[inode]=99; // <--
    for (ipar=0;ipar<node_parm[inode].size();++ipar){
      int parm=node_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      node_quan[inode].push_back(quan); // <--
      if (quan<first_quan[inode]) first_quan[inode]=quan; // <--
    }
  }
}
// ----------------------------------------------------------------------
void get_num_den_part(int nnode, vector<int> * node_parm,
		      int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		      int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef,
		      vector<int> baseparm, double* baseseed,
		      double* node_num, // output
		      double* node_den, // output
		      double* node_val, // output
		      vector<double> * node_part// output
		      ){
  int ipar;
  //
  // General Case:
  // R(P_i) = (sum_j=1,n alpha_j P_j) / (sum_k=1,m beta_k P_k)
  // dR/dP_i= (1/sum_k=1,m beta_i P_i) * alpha_j * delta_ij - (R/sum_k=1,m beta_k P_k) * beta_k * delta_ik
  //
  // Special Case 1:
  // R(P_i) = (sum_j=1,n alpha_j P_j) 
  // dR/dP_i= alpha_i 
  //
  // Special Case 2: 
  // R(P_i) = (sum_j=1,n alpha_j P_j) / (sum_k=1,m beta_k P_k)
  // Case 2a: if P_i only in numerator: dR/dP_i = (1/sum_k=1,m beta_k P_k) * alpha_i
  // Case 2b: if P_i only in denominator: dR/dP_i = - (R/sum_k=1,m beta_k P_k) * beta_i
  // Case 2c: if P_i both in numerator and denominator: dR/dP_i = (1/sum_k=1,m beta_i P_i) * alpha_i - (R/sum_k=1,m beta_k P_k) * beta_i
  //
  for (int inode=0;inode<nnode;++inode){
    //
    double numerator=0; // sum_j=1,n alpha_j P_j
    if (node_num_npar[inode]>0) {
      for (ipar=0;ipar<node_num_npar[inode];++ipar){
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	numerator+=(node_num_coef[inode].at(ipar))*(baseseed[quan-1]); 
      }
    }
    node_num[inode]=numerator; // <--
    //
    double denominator=0; // sum_k=1,m beta_k P_k
    if (node_den_npar[inode]>0) {
      for (ipar=0;ipar<node_den_npar[inode];++ipar){
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	denominator+=(node_den_coef[inode].at(ipar))*(baseseed[quan-1]); 
      }
    }
    node_den[inode]=denominator; // <--
    //
    node_val[inode]=(denominator!=0) ? numerator/denominator : numerator ; // <--
    //
    for (ipar=0;ipar<node_parm[inode].size();++ipar){
      int parm=node_parm[inode].at(ipar);
      //
      vector<int>::iterator it_num=find(node_num_parm[inode].begin(),node_num_parm[inode].end(),parm); bool is_in_num = it_num != node_num_parm[inode].end();
      vector<int>::iterator it_den=find(node_den_parm[inode].begin(),node_den_parm[inode].end(),parm); bool is_in_den = it_den != node_den_parm[inode].end();
      double partial=0;
      if (node_den_npar[inode]==0 && is_in_num) {//case 1
	partial = node_num_coef[inode].at(it_num - node_num_parm[inode].begin());
      } else if ( is_in_num && !is_in_den) { // case 2a
	partial = (1./denominator) * (node_num_coef[inode].at(it_num - node_num_parm[inode].begin()));
      } else if (!is_in_num &&  is_in_den) { // case 2b
	partial = -1. * (numerator/(denominator*denominator)) * (node_den_coef[inode].at(it_den - node_den_parm[inode].begin()));
      } else if ( is_in_num &&  is_in_den) { // case 2c
	partial = (1./denominator) * (node_num_coef[inode].at(it_num - node_num_parm[inode].begin())) -1. * (numerator/(denominator*denominator)) * (node_den_coef[inode].at(it_den - node_den_parm[inode].begin()));
      }
      node_part[inode].push_back(partial); // <--
      //      if (inode==14) cout << " inode = " << inode << " ipar = " << ipar << " parm = " << parm 
      //			  << " numerator = " << numerator << " denominator = " << denominator << " partial = " << partial 
      //      			  << " it_num - node_num_parm[inode].begin() = " << it_num - node_num_parm[inode].begin() 
      //      			  << endl;
    }
  }
}
// ----------------------------------------------------------------------
void combine(int uconstrain, 
	     int nmeas, int* weak, int weakcompare, int& nmeas_noweak,
	     int* measnode, double* measvalue, TMatrixD InvMeasErrorMatrix,
	     int* node_num_npar, int* node_den_npar, vector<int> * node_quan, vector<int> * node_parm, bool* node_is_base,
	     int nbase, vector<int> basegamma,  
	     double* baseseed, double* node_num, double* node_den, vector<double> * node_part,
	     double* basevalue_fit, // output
	     double* baseerror_fit, // output
	     double** basecov_fit,  // output
	     double** basecorr_fit, // output
	     double&_chisquared     // output
	     ){
  //
  // measuredquantity = f(wxyz) = f(0) + f'w (w-w0) + f'x (x-x0) + f'y (y-y0) + f'z (z-z0)
  //                            = [f(0) - f'w w0 - f'x x0 - f'y y0 - f'z z0] + f'w w + f'x x + f'y y + f'z z

  // [measuredvalue - measuredquantity] = [measuredvalue - (f(0) - f'w w0 - f'x x0 - f'y y0 - f'z z0) - (f'w w + f'x x + f'y y + f'z z) ] 
  //                                    = Xvector + Delta^T Vvector
  // Xvector = vector [dimension = nmeas] of adjusted measured value. In this case, adjust measurement = M to M - f(0) + f'w w0 + f'x x0 + f'y y0 + f'z z0
  // Delta^T = matrix of dimension (nmeas,nbase). In this case it has dimension 1x4 = (-f'w -f'x -f'y -f'z)
  // Vvector = vector [dimension = nbase] of base quantity. In this case Vvector^T is 1x4 = (w,x,y,z)
  //
  int i,j,imeas,jmeas,inode,ibase,jbase,ipar;
  //
  nmeas_noweak=0; 
  for (imeas=0;imeas<nmeas;++imeas) {
    if (weak[imeas]==weakcompare) continue; 
    ++nmeas_noweak;
  }
  //
  TMatrixD InvMeasErrorMatrix_noweak(nmeas_noweak,nmeas_noweak);
  i=-1;
  for (imeas=0;imeas<nmeas;++imeas) {
    if (weak[imeas]==weakcompare) continue; 
    ++i;
    j=-1;
    for (jmeas=0;jmeas<nmeas;++jmeas) {
      if (weak[jmeas]==weakcompare) continue; 
      ++j;
      InvMeasErrorMatrix_noweak[i][j]=InvMeasErrorMatrix[imeas][jmeas];
    }
  }
  //
  int nbase_u = (uconstrain) ? nbase-1 : nbase;
  TMatrixD Xvector(nmeas_noweak,1);
  TMatrixD Delta(nbase_u,nmeas_noweak);
  TMatrixD DeltaT(nmeas_noweak,nbase_u);
  TMatrixD SDInvMatrix(nbase_u,nbase_u);
  TMatrixD SDMatrix(nbase_u,nbase_u);
  Double_t sd_det;
  TMatrixD Vvector(nbase_u,1);
  TMatrixD XPlusDTV(nmeas_noweak,1);
  TMatrixD XPlusDTVT(1,nmeas_noweak);
  TMatrixD ChiSquared(1,1);
  //
  i=-1;
  for (imeas=0;imeas<nmeas;++imeas){
    if (weak[imeas]==weakcompare) continue;
    ++i;
    //
    for (ibase=0;ibase<nbase_u;++ibase) Delta[ibase][i] = 0;
    //
    inode=measnode[imeas];
    if (uconstrain &&  // SPECIAL CASE [because these nodes contain Gamma103 [used to express unitarity constraint] [N.B. Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)] 
	(inode==N_GAMMA102 || inode==N_GAMMA103)) { 
      Xvector[i][0] = measvalue[imeas] - 1 ;
      for (ibase=0;ibase<nbase_u;++ibase) {
	Delta[ibase][i] = 1;
	if (inode==N_GAMMA102 && ibase==M_GAMMA104) Delta[ibase][i] = 0; // remove dependence on ibase=Gamma104 for inode=Gamma102
      }
    } else {
      Xvector[i][0] = measvalue[imeas];
      double offset = -node_num[inode]; if (node_den_npar[inode]>0) offset /= node_den[inode];
      if (!node_is_base[inode]) {
	Xvector[i][0]+=offset;
      }
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int quan=node_quan[inode].at(ipar);
	double partial=node_part[inode].at(ipar);
	if (!node_is_base[inode]) {
	  Xvector[i][0] += partial*baseseed[quan-1]; 
	}
	Delta[quan-1][i] = -partial;
      }
    }
  }
  //
  DeltaT.Transpose(Delta);
  SDInvMatrix = Delta * InvMeasErrorMatrix_noweak * DeltaT;
  SDMatrix = SDInvMatrix;
  SDMatrix.Invert(&sd_det);
  Vvector = SDMatrix * Delta * InvMeasErrorMatrix_noweak * Xvector; Vvector*=-1;
  XPlusDTV = Xvector + DeltaT*Vvector;
  XPlusDTVT.Transpose(XPlusDTV);
  ChiSquared = XPlusDTVT * InvMeasErrorMatrix_noweak * XPlusDTV;
  _chisquared = ChiSquared[0][0];
  //
  for (ibase=0;ibase<nbase_u;++ibase) {
    basevalue_fit[ibase] = Vvector[ibase][0] ;
    baseerror_fit[ibase] =  TMath::Sqrt(SDMatrix[ibase][ibase]);
  }
  if (uconstrain) {
    basevalue_fit[nbase-1] = 1;
    double error_squared = 0;
    for (ibase=0;ibase<nbase_u;++ibase) {
      basevalue_fit[nbase-1] -= Vvector[ibase][0] ;
      for (jbase=0;jbase<nbase_u;++jbase) {
	error_squared += SDMatrix[ibase][jbase];
      }
      baseerror_fit[nbase-1] = TMath::Sqrt(error_squared);
    }
  }
  //
  for (ibase=0;ibase<nbase_u;++ibase) {
    for (jbase=0;jbase<nbase_u;++jbase) {
      basecov_fit [ibase][jbase] = SDMatrix[ibase][jbase];
      basecorr_fit[ibase][jbase] = (ibase==jbase) ? 1 : basecov_fit[ibase][jbase]/(baseerror_fit[ibase]*baseerror_fit[jbase]);
    }
  }
  if (uconstrain) {
    for (ibase=0;ibase<nbase_u;++ibase){
      basecov_fit[ibase][nbase-1] = 0;
      basecov_fit[nbase-1][ibase] = 0;
      for (jbase=0;jbase<nbase_u;++jbase) {
	basecov_fit[ibase][nbase-1] -= basecov_fit[ibase][jbase];
	basecov_fit[nbase-1][ibase] -= basecov_fit[jbase][ibase];
      }
    }
    basecov_fit [nbase-1][nbase-1] = baseerror_fit[nbase-1]*baseerror_fit[nbase-1];
    basecorr_fit[nbase-1][nbase-1] = 1;
    for (ibase=0;ibase<nbase_u;++ibase){
      basecorr_fit[ibase][nbase-1] = basecov_fit[ibase][nbase-1] / (baseerror_fit[ibase]*baseerror_fit[nbase-1]);
      basecorr_fit[nbase-1][ibase] = basecov_fit[ibase][nbase-1] / (baseerror_fit[ibase]*baseerror_fit[nbase-1]);
    }
  }
}
// ----------------------------------------------------------------------
void process_measurements(int nmeas, int* measnode, double* measerror, double** corrmat,
			  TMatrixD & MeasErrorMatrix,  // output
			  TMatrixD & InvMeasErrorMatrix, //output
			  vector<int> *icorrj, //output
			  int &ncorrij, //output
			  int *ifirstj, //output
			  vector<int> *veccorrij, //output
			  int & newnode_all, //output
			  int* nodegroup_all, //output
			  int* ncycle, // output
			  int* n_nodegroup, // output
			  int* ngroup //output
			  ){
  //
  int i,j,inode;
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) {
	MeasErrorMatrix(i,j)=measerror[i]*measerror[i];
      } else {
	MeasErrorMatrix(i,j)=corrmat[i][j]*measerror[i]*measerror[j]; // by construction diagonal elements of corrmat has not been filled
      }
    }
  }
  //
  InvMeasErrorMatrix = MeasErrorMatrix;
  Double_t determinant;
  InvMeasErrorMatrix.Invert(&determinant);  
  //  cout << Form("Determinant of MeasErrorMatrix = %10.5g\n", determinant);
  //
  for (i=0;i<nmeas;++i){
    for (j=0;j<nmeas;++j){
      if (InvMeasErrorMatrix[i][j]!=0.0) {
	icorrj[i].push_back(j);
      }
    }
  }
  //
  ncorrij=0;
  for (i=0;i<nmeas;++i) {
    if (icorrj[i].size() > 1) {
      bool isnew=true;
      for (j=0;j<i;++j){
	if (icorrj[j].size() > 1) {
	  if (find(icorrj[j].begin(),icorrj[j].end(),i)!=icorrj[j].end()) {
	    isnew=false;
	    break;
	  }
	}
      }
      if (isnew) {
	ifirstj[ncorrij] = i;
	++ncorrij;
      }
      //      cout << "i+1,isnew,ncorrij = " << i+1 << " " << isnew << " " << ncorrij << endl;
    }
  }
  //
  for (i=0;i<ncorrij;++i){
    veccorrij[i].insert(veccorrij[i].end(),icorrj[ifirstj[i]].begin(),icorrj[ifirstj[i]].end());
  }
  //
  newnode_all=-1;
  vector<int> vector_measnodes_all;
  for (i=0;i<nmeas;++i) {
    if (icorrj[i].size()>1) continue;
    ncycle[i] = -1; // cycle number for uncorrelated measurements
    inode=measnode[i];
    vector<int>::iterator itt=find(vector_measnodes_all.begin(),vector_measnodes_all.end(),inode);
    bool is_newnode_all = itt == vector_measnodes_all.end();
    if (is_newnode_all) {
      ++newnode_all;
      nodegroup_all[i] = newnode_all; // node-group number for new uncorrelated measurements
    } else {
      for (j=0;j<i;++j) {
	if (icorrj[j].size()>1) continue;
	if (measnode[j]==inode) {
	  nodegroup_all[i] = nodegroup_all[j]; // node-group number for old uncorrelated measurements
	  break;
	} 
      }
    }
    vector_measnodes_all.push_back(inode);
    //    cout << " i+1 = " << i+1 << " inode = " << inode << " is_newnode_all = " << is_newnode_all << " newnode_all = " << newnode_all << " nodegroup = " << nodegroup_all[i] << endl;
  }
  for (i=0;i<ncorrij;++i) {
    ++newnode_all;
    for (j=0;j<veccorrij[i].size();++j) {
      int meas=veccorrij[i][j];
      ncycle[meas] = i; // cycle number for each correlated measurement
      nodegroup_all[meas] = newnode_all; // node-group number for correlated measurements
    }
  }
  for (inode=0;inode<=newnode_all;++inode){
    n_nodegroup[inode]=0;
    for (i=0;i<nmeas;++i){
      if (nodegroup_all[i]==inode){
	n_nodegroup[inode]+=1;
      }
    }
  }
  for (inode=0;inode<=newnode_all;++inode){
    for (i=0;i<nmeas;++i){
      if (nodegroup_all[i]==inode) {
	ngroup[i] = n_nodegroup[nodegroup_all[i]];
      }
    }
  }
}
// ----------------------------------------------------------------------
void print_measinfo(FILE* thisfile,
		    int nmeas, int* measnode, double* measvalue, double* measerror, double** corrmat,
		    string* expname, string* meastitle, string* measgammaname,
		    char** a_nodename, int* node_num_npar, int* node_den_npar, vector<int> * node_quan, bool* node_is_base,
		    int nbase, vector<int> baseparm, vector<int> basegamma, string* basenodename, string* basetitle,
		    int newnode_all, int* nodegroup_all, int* ncycle, int ncorrij, int* ifirstj, vector<int> *icorrj, vector<int> *veccorrij, int* ngroup){
  //
  int i, inode, ibase, j;
  //
  fprintf (thisfile, "nmeas = %d \n", nmeas);
  //
  vector<int> basequan_used_in_measured_basenodes;
  vector<int> basequan_used_in_measured_derivednodes;
  for (i=0;i<nmeas;++i){
    inode=measnode[i];
    if (!node_is_base[inode]) {
      basequan_used_in_measured_derivednodes.insert(basequan_used_in_measured_derivednodes.end(),node_quan[inode].begin(),node_quan[inode].end());
      sort(basequan_used_in_measured_derivednodes.begin(),basequan_used_in_measured_derivednodes.end());
      vector<int>::iterator new_end=unique(basequan_used_in_measured_derivednodes.begin(),basequan_used_in_measured_derivednodes.end());
      basequan_used_in_measured_derivednodes.erase(new_end,basequan_used_in_measured_derivednodes.end());
    } else { // base node
      basequan_used_in_measured_basenodes.insert(basequan_used_in_measured_basenodes.end(),node_quan[inode].begin(),node_quan[inode].end());
      sort(basequan_used_in_measured_basenodes.begin(),basequan_used_in_measured_basenodes.end());
      vector<int>::iterator new_end=unique(basequan_used_in_measured_basenodes.begin(),basequan_used_in_measured_basenodes.end());
      basequan_used_in_measured_basenodes.erase(new_end,basequan_used_in_measured_basenodes.end());
    }
  }
  fprintf (thisfile, "basequan_used_in_measured_basenodes.size() = %d ", basequan_used_in_measured_basenodes.size());
  fprintf (thisfile, "[IF THIS IS NOT = %d, IT MEANS SOME BASE QUANTITIES ARE NOT MEASURED DIRECTLY BUT ONLY VIA MEASUREMENTS OF DERIVED NODES]\n", nbase); 
  for (ibase=0;ibase<nbase;++ibase) {
    vector<int>::iterator it=find(basequan_used_in_measured_basenodes.begin(),basequan_used_in_measured_basenodes.end(),ibase+1); 
    bool not_present = it == basequan_used_in_measured_basenodes.end();
    if (not_present) fprintf (thisfile, "BASE = %d GAMMA = %d PARM = %d NODE = %s TITLE = %s not measured directly \n",
			      ibase+1, basegamma[ibase], baseparm[ibase], basenodename[ibase].data(), basetitle[ibase].data());
  }
  fprintf (thisfile, "\n");
  fprintf (thisfile, "basequan_used_in_measured_derivednodes.size() = %d \n", basequan_used_in_measured_derivednodes.size());
  fprintf (thisfile, "\n");
  //  for (ibase=0;ibase<basequan_used_in_measured_derivednodes.size();++ibase) {
  //    int quan = basequan_used_in_measured_derivednodes.at(ibase) ;
  for (int quan=1;quan<nbase+1;++quan) {
    fprintf (thisfile, "basetitle = %s with quan = %d gamma = %d appears in ",basetitle[quan-1].data(), quan, basegamma[quan-1]);
    int iquan_occurance=0;
    int inode_occurance=0;
    vector<int> vector_measnodes;
    for (i=0;i<nmeas;++i){
      inode=measnode[i];
      vector<int>::iterator itt=find(vector_measnodes.begin(),vector_measnodes.end(),inode);
      bool is_newnode = itt == vector_measnodes.end();
      for (int ii=0;ii<node_quan[inode].size();++ii){
	if (quan == node_quan[inode].at(ii)) { 
	  ++iquan_occurance;
	  if (is_newnode) {
	    ++inode_occurance;
	  }
	}
      }
      vector_measnodes.push_back(inode);
    }
    fprintf (thisfile, "%d measurements of %d nodes. \nThese %d measurements are :\n", iquan_occurance, inode_occurance, iquan_occurance);
    iquan_occurance = 0;
    vector<int> vector_of_measurements_correlated_with_other_measurements;
    for (i=0;i<nmeas;++i){
      bool meas_with_correlation = false;
      inode=measnode[i];
      for (int ii=0;ii<node_quan[inode].size();++ii){
	if (quan == node_quan[inode].at(ii)) { 
	  ++iquan_occurance;
	  fprintf (thisfile, " (%d) IMEAS = %d NODE = %d NAME = %s GAMMA = %s TITLE = %s EXP = %s\n",
		   iquan_occurance, i+1, inode, a_nodename[inode], measgammaname[i].data(), meastitle[i].data(), expname[i].data());
	  for (j=0;j<nmeas;++j){
	    if (corrmat[i][j] != 0) {
	      meas_with_correlation = true;
	      break;
	    }
	  }
	}
      }
      if (meas_with_correlation) vector_of_measurements_correlated_with_other_measurements.push_back(i+1); // <- NOTE: +1
    }
    if (vector_of_measurements_correlated_with_other_measurements.size()!=0) {
      fprintf (thisfile, "vector_of_measurements_correlated_with_other_measurements has size = %d : ",vector_of_measurements_correlated_with_other_measurements.size());
      for (int i=0; i < vector_of_measurements_correlated_with_other_measurements.size(); ++i) {
	int imeas = vector_of_measurements_correlated_with_other_measurements[i] -1; // <- NOTE: -1
	fprintf (thisfile, "%d(%s) ",imeas+1, expname[imeas].data()); 
      }
      fprintf (thisfile, "\n");
    } else {
      fprintf (thisfile, "vector_of_measurements_correlated_with_other_measurements has size = 0\n");
    }
    fprintf (thisfile, "\n");
  }
  //
  fprintf (thisfile, "ncorrij = %d \n",ncorrij);
  for (i=0;i<ncorrij;++i){
    fprintf (thisfile, "i = %d ifirstj[i]+1 = %3d ",i,ifirstj[i]+1);
    fprintf (thisfile, "veccorrij[i].size() = %2d veccorrij[i][j]+1 = ",veccorrij[i].size());
    for (j=0;j<veccorrij[i].size();++j) fprintf (thisfile, "%d ",veccorrij[i][j]+1);
    fprintf (thisfile, "GammaName = "); 
    for (j=0;j<veccorrij[i].size();++j) fprintf (thisfile, "%s ",measgammaname[veccorrij[i][j]].data());
    fprintf (thisfile,"\n");
  }
  //
  for (inode=0;inode<=newnode_all;++inode){
    for (i=0;i<nmeas;++i){
      if (nodegroup_all[i]==inode) {
	fprintf (thisfile, "i+1 = %3d nodegroup = %2d ngroup = %2d icorrj = %2d ncycle = %2d meas = %10.5e +- %10.5e %6s %s \n",
		 i+1,nodegroup_all[i],ngroup[i],icorrj[i].size(), ncycle[i], measvalue[i], measerror[i], expname[i].data(), meastitle[i].data());
      }
    }
  }
}
// ----------------------------------------------------------------------
void print_measfile(FILE* thisfile, int p, int uconstrain,
		    int nmeas, int* weak, int weakcompare, string* measgammaname, int* measnode, 
		    string* expname, string* author, string* year, string* meastitle,
		    double* measvalue, double* measerror, double** corrmat,
		    char** a_nodename, vector<int> * node_parm,
		    int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		    int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef, 
		    int nbase, vector<int> baseparm, vector<int> basegamma, string* basetitle, int* first_quan, double* baseseed){
  int inode,ipar;
  int iimeas=0;
  for (int i=0;i<nmeas;++i){
    if (weak[i]==weakcompare) continue;
    ++iimeas;
    //
    fprintf (thisfile, "* IMEAS = %d \n",     iimeas);
    fprintf (thisfile, "* GAMMANAME = %s \n", measgammaname[i].data());
    fprintf (thisfile, "* DECAYNAME = %s \n", meastitle[i].data());
    //
    inode=measnode[i];
    //    cout << "i = " << i << " inode = " << inode << " nodename = " << a_nodename[inode]
    //	 << " node_num_npar[inode] = " << node_num_npar[inode] << " node_num_parm[inode].size() = " << node_num_parm[inode].size()  
    //	 << " node_den_npar[inode] = " << node_den_npar[inode] << " node_den_parm[inode].size() = " << node_den_parm[inode].size()  
    //	 << " node_parm[inode].size() = " << node_parm[inode].size() 
    //	 << endl;
    //
    fprintf (thisfile, 
	     "* NODENAME = %s found at inode+1 = %d has %d, %d, %d base quantities in numerator, denominator and both [excluding overlap] \n", 
	     a_nodename[inode], inode+1, node_num_npar[inode], node_den_npar[inode],node_parm[inode].size());
    //
    for (ipar=0;ipar<node_num_parm[inode].size();++ipar){
      int parm=node_num_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      fprintf (thisfile,"*                numerator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %20.10f coef = %20.10f\n",
	       inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_num_coef[inode].at(ipar));
    }
    for (ipar=0;ipar<node_den_parm[inode].size();++ipar){
      int parm=node_den_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      fprintf (thisfile,"*              denominator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %20.10f coef = %20.10f\n",
	       inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_den_coef[inode].at(ipar));
    }
    fprintf (thisfile, "*  first quantity measured by inode+1 = %d has gamma = %d parm = %d quan = %d title = %s\n",
	     inode+1,basegamma[first_quan[inode]-1],baseparm[first_quan[inode]-1],first_quan[inode],basetitle[first_quan[inode]-1].data());
    //
    fprintf (thisfile, "\nBEGIN %s Gamma%s pub.%s.%s \n\n", expname[i].data(), measgammaname[i].data(), author[i].data(), year[i].data());
    if (p==0) {//COMBOS
      if (uconstrain &&  // SPECIAL CASE [because these nodes contain Gamma103 [used to express unitarity constraint] [N.B. Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)] 
	  (inode==N_GAMMA102 || inode==N_GAMMA103)) { 
	fprintf (thisfile, "MEASUREMENT  m_Gamma%d statistical systematic \n",3);
	fprintf (thisfile, "DATA         m_Gamma%d statistical systematic \n",3);
      } else {
	fprintf (thisfile, "MEASUREMENT  m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	fprintf (thisfile, "DATA         m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
      }
    } else if (p==1) {//ALUCOMB
      fprintf (thisfile, "MEASUREMENT  m_Gamma%s statistical systematic \n",measgammaname[i].data());
      fprintf (thisfile, "DATA         m_Gamma%s statistical systematic \n",measgammaname[i].data());
    }
    fprintf (thisfile, "             %20.10f %20.10f  0 \n",measvalue[i],measerror[i]);
    bool firstcorr=true;
    int jjmeas=0;
    for (int j=0;j<nmeas;++j) {
      if (weak[j]==weakcompare) continue;
      ++jjmeas;
      if (corrmat[i][j]!=0) {
	if (firstcorr) {fprintf (thisfile, " \n"); firstcorr=false;}
	fprintf (thisfile, "STAT_CORR_WITH %s Gamma%s pub.%s.%s %20.10f ! IMEAS = %d \n", 
		 expname[j].data(), measgammaname[j].data(), author[j].data(), year[j].data(), corrmat[i][j], jjmeas);
      }
    }
    fprintf (thisfile, " \nEND \n\n");
  }
}
// ----------------------------------------------------------------------
void print_avefile(FILE* thisfile, int p, int uconstrain,
		   int nmeas, int* weak, int weakcompare, string* measgammaname, int* measnode,
		   vector<int> * node_parm, vector<int> * node_quan, bool* node_is_base,
		   int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		   int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef, 
		   int nbase, vector<int> baseparm, vector<int> basegamma, string* basetitle, int* first_quan, 
		   double* baseseed, double* node_num, double* node_den,  vector<double> * node_part){
  int i,inode,ipar,iimeas,ibase;
#if defined USING_NBASE31 || defined USING_NBASE39
  fprintf (thisfile, "BEGIN   PDG-BABAR-BELLE all_methods \n\n");
#elif defined USING_NBASE40
  fprintf (thisfile, "BEGIN   PDG+BABAR+BELLE all_methods \n\n");
#endif
  fprintf (thisfile, "COMBINE * * * \n\n");
  for (ibase=0;ibase<nbase;++ibase){
    if (p==0&&uconstrain&&ibase==(nbase-1)){/* skip */} else {
      fprintf (thisfile, "MEASUREMENT m_Gamma%d statistical systematic   ! NQUAN = %d \n",basegamma.at(ibase),ibase+1);  
    }
  }
  if (p==0){
    fprintf (thisfile, "\nCALL DUMP_MASTER_INC \n\n");
    fprintf (thisfile, "SPARAMETER CHI2_N_SYM_PRT -1.0 0. \n\n");
    fprintf (thisfile, "SPARAMETER CHI2_N_SYM_INV 0 0    \n\n");
  }
  //
  vector<int> vector_measnodes;
  int usum=0;//number of [unique] measurements to be expressed as linearized sum of base quantities
  for (i=0;i<nmeas;++i){
    if (weak[i]==weakcompare) continue;
    inode=measnode[i];
    vector<int>::iterator itt=find(vector_measnodes.begin(),vector_measnodes.end(),inode);
    bool is_newnode = itt == vector_measnodes.end();
    if (!node_is_base[inode]) {
      if (p==1&&is_newnode){//new node
	++usum;
	fprintf (thisfile, "MEASUREMENT m_Gamma%s statistical systematic   ! NQUAN = %d \n",measgammaname[i].data(),ibase+usum);  
      }
    }
    vector_measnodes.push_back(inode);      
  }
  //
  if (p==1&&!uconstrain){
    fprintf (thisfile, "*--- dummy node to remove unitarity constraint\nMEASUREMENT m_Gamma998 statistical systematic\n");
  }
  if (p==1) {
    fprintf (thisfile, "*--- Gamma110, only present in constraint\nMEASUREMENT m_Gamma110 statistical systematic\n");
  }
  //
  int isum=0;//number of          measurements to be expressed as linearized sum of base quantities
  vector_measnodes.erase(vector_measnodes.begin(),vector_measnodes.end());
  iimeas=0;
  for (int i=0;i<nmeas;++i){
    if (weak[i]==weakcompare) continue;
    ++iimeas;
    inode=measnode[i];
    vector<int>::iterator itt=find(vector_measnodes.begin(),vector_measnodes.end(),inode);
    bool is_newnode = itt == vector_measnodes.end();
    if (p==1 && !is_newnode) continue; // ALUCOMB needs it only once
    if (!node_is_base[inode]) {
      //
      if (p==0 && (inode==N_GAMMA102 || inode==N_GAMMA103)) continue; // SPECIAL CASE [because these are derived nodes containing Gamma103]
      ++isum; // translate C index to Fortran index
      //
      // PRINT NODE DEFINITION
      //
      fprintf (thisfile, "\n* Gamma%s = ",measgammaname[i].data());
      for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
	if (ipar==0) { fprintf (thisfile, "(") ; } else {fprintf (thisfile, " + ") ;}
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (thisfile, "%20.10f*Gamma%d",node_num_coef[inode].at(ipar), basegamma[quan-1]);
	if (ipar==node_num_parm[inode].size()-1) fprintf (thisfile, ")");
      }
      if (node_den_parm[inode].size()==0) fprintf (thisfile, "\n") ; 
      for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
	if (ipar==0) { fprintf (thisfile, " / (") ; } else {fprintf (thisfile, " + ") ;}
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (thisfile, "%20.10f*Gamma%d",node_den_coef[inode].at(ipar), basegamma[quan-1]);
	if (ipar==node_den_parm[inode].size()-1) fprintf (thisfile, ")\n");
      }
      //
      double offset = -node_num[inode]; // - [ f(x0,y0) - df/dx(x=x0) x0 - df/dy(y=y0) y0 - ...]
      if (node_den_npar[inode]>0) offset /= node_den[inode];
      //	  cout << inode << " " << node_num[inode] << " " << node_den[inode] << " " << offset << " " << measvalue[i] << endl;
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int quan=node_quan[inode].at(ipar);
	double partial=node_part[inode].at(ipar);
	offset += partial*baseseed[quan-1]; 
      }
      if (p==0) { // COMBOS
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %d\n",isum,iimeas,node_parm[inode].size()); 
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD %20.10f 1.0\n",isum,offset); 
      } else if (p==1) { // ALUCOMB
	fprintf (thisfile, "CONSTRAINT Gamma%s.c %20.10f Gamma%s -1", measgammaname[i].data(), offset, measgammaname[i].data());
      }
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int quan=node_quan[inode].at(ipar);
	double partial=node_part[inode].at(ipar);
	if (p==0) { // COMBOS
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %20.10f ! Gamma%d \n",isum,ipar+1,quan,partial,basegamma.at(quan-1));
	} else if (p==1) { // ALUCOMB
	  fprintf (thisfile, " Gamma%d %20.10f", basegamma[quan-1], partial);
	}
      }
      if (p==1) fprintf (thisfile, "\n");
    }
    vector_measnodes.push_back(inode);      
  }
  if (p==0) {
    iimeas=0;
    for (int i=0;i<nmeas;++i){
      if (weak[i]==weakcompare) continue;
      ++iimeas;
      inode=measnode[i];
      //
      if (inode==N_GAMMA102) { // SPECIAL CASE [because these are derived nodes containing Gamma103]
	if (!uconstrain) {
	  fprintf (thisfile, "\n*Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)\n");
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %d\n",++isum,iimeas,node_parm[inode].size()); 
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD %20.10f 1.0\n",isum,0); 
	  for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	    int quan=node_quan[inode].at(ipar);
	    double partial=node_part[inode].at(ipar);
	    fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %20.10f ! Gamma%d \n",isum,ipar+1,quan,partial,basegamma.at(quan-1));
	  }
	} else {
	  fprintf (thisfile, "\n*Gamma102 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	  fprintf (thisfile, "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	  fprintf (thisfile, "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
#if defined USING_NBASE31
	  fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	  fprintf (thisfile, "*             - Gamma94 - Gamma126 - Gamma128 - Gamma150 - Gamma152\n");
#elif defined USING_NBASE39
	  fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma802 - Gamma803 - Gamma93\n");
	  fprintf (thisfile, "*             - Gamma94 - Gamma126 - Gamma128 - Gamma800 - Gamma151 - Gamma152\n");
	  fprintf (thisfile, "*             - Gamma130 - Gamma132 - Gamma44 - Gamma53\n");
	  fprintf (thisfile, "*             - Gamma49  - Gamma804 - Gamma805\n");
#elif defined USING_NBASE40
	  fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma802 - Gamma803 - Gamma93\n");
	  fprintf (thisfile, "*             - Gamma94 - Gamma126 - Gamma128 - Gamma800 - Gamma151 - Gamma152\n");
	  fprintf (thisfile, "*             - Gamma130 - Gamma132 - Gamma44 - Gamma53 - Gamma801\n");
	  fprintf (thisfile, "*             - Gamma49  - Gamma804 - Gamma805\n");
#endif
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %2d \n",++isum,iimeas,nbase-2); 
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma102; thats why the coefficients below have - sign 
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
#if defined USING_NBASE31
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
#elif defined USING_NBASE39 || defined USING_NBASE40
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma802 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma803 \n",isum);
#endif
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_26 27 -1 ! Gamma126\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_27 28 -1 ! Gamma128\n",isum);
#if defined USING_NBASE31
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_28 29 -1 ! Gamma150\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_29 30 -1 ! Gamma152\n",isum);
#elif defined USING_NBASE39 || defined USING_NBASE40
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_28 29 -1 ! Gamma800\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_29 30 -1 ! Gamma151\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_30 31 -1 ! Gamma152\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_31 32 -1 ! Gamma130\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_32 33 -1 ! Gamma132\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_33 34 -1 ! Gamma44\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_34 35 -1 ! Gamma53\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_35 36 -1 ! Gamma49\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_36 37 -1 ! Gamma804\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_37 38 -1 ! Gamma805\n",isum);
#if defined USING_NBASE40
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_38 39 -1 ! Gamma801\n",isum);
#endif
#endif
	}
      }
      //
      if (uconstrain && inode==N_GAMMA103) { // SPECIAL CASE [because Gamma103 is used to express unitarity constraint]
	fprintf (thisfile, "\n*Gamma103 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	fprintf (thisfile, "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	fprintf (thisfile, "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
#if defined USING_NBASE31
	fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	fprintf (thisfile, "*             - Gamma94 - Gamma104 - Gamma126 - Gamma128 - Gamma150 - Gamma152\n");
#elif defined USING_NBASE39
	fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma802 - Gamma803 - Gamma93\n");
	fprintf (thisfile, "*             - Gamma94 - Gamma104 - Gamma126 - Gamma128 - Gamma800 - Gamma151 - Gamma152\n");
	fprintf (thisfile, "*             - Gamma130 - Gamma132 - Gamma44 - Gamma53\n");
	fprintf (thisfile, "*             - Gamma49  - Gamma804 - Gamma805\n");
#elif defined USING_NBASE40
	fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma802 - Gamma803 - Gamma93\n");
	fprintf (thisfile, "*             - Gamma94 - Gamma104 - Gamma126 - Gamma128 - Gamma800 - Gamma151 - Gamma152\n");
	fprintf (thisfile, "*             - Gamma130 - Gamma132 - Gamma44 - Gamma53 - Gamma801\n");
	fprintf (thisfile, "*             - Gamma49  - Gamma804 - Gamma805\n");
#endif
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %2d \n",++isum,iimeas,nbase-1); 
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma103; thats why the coefficients below have - sign 
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
#if defined USING_NBASE31
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
#elif defined USING_NBASE39 || defined USING_NBASE40
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma802 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma803 \n",isum);
#endif
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_26 26 -1 ! Gamma104\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_27 27 -1 ! Gamma126\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_28 28 -1 ! Gamma128\n",isum);
#if defined USING_NBASE31
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_29 29 -1 ! Gamma150\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_30 30 -1 ! Gamma152\n",isum);
#elif defined USING_NBASE39 || defined USING_NBASE40
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_29 29 -1 ! Gamma800\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_30 30 -1 ! Gamma151\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_31 31 -1 ! Gamma152\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_32 32 -1 ! Gamma130\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_33 33 -1 ! Gamma132\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_34 34 -1 ! Gamma44\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_35 35 -1 ! Gamma53\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_36 36 -1 ! Gamma49\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_37 37 -1 ! Gamma804\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_38 38 -1 ! Gamma805\n",isum);
#if defined USING_NBASE40
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_39 39 -1 ! Gamma801\n",isum);
#endif
#endif
      }
    }
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_NSUM  %d 0 \n",isum); 
  }
  if (p==1) {
    if (uconstrain) {
      fprintf (thisfile, "\n* unitarity constraint applied (sum of base nodes without dummy node)\n");
    } else {
      fprintf (thisfile, "\n* unitarity constraint NOT applied (sum of base nodes with dummy node)\n");
    }
    fprintf (thisfile, "CONSTRAINT GammaAll 1\n");
    fprintf (thisfile, "  Gamma3   1 Gamma5   1 Gamma9   1 Gamma10  1 Gamma14  1 Gamma16  1\n");
    fprintf (thisfile, "  Gamma20  1 Gamma23  1 Gamma27  1 Gamma28  1 Gamma30  1 Gamma35  1\n");
    fprintf (thisfile, "  Gamma37  1 Gamma40  1 Gamma42  1 Gamma47  1 Gamma48  1 Gamma62  1\n");
#if defined USING_NBASE31
    fprintf (thisfile, "  Gamma70  1 Gamma77  1 Gamma78  1 Gamma85  1 Gamma89  1 Gamma93  1\n");
    fprintf (thisfile, "  Gamma94  1 Gamma103 1 Gamma104 1 Gamma126 1 Gamma128 1 Gamma150 1 Gamma152 1\n");
#elif defined USING_NBASE39
    fprintf (thisfile, "  Gamma70  1 Gamma77  1 Gamma78  1 Gamma802 1 Gamma803   Gamma93  1\n");
    fprintf (thisfile, "  Gamma94  1 Gamma103 1 Gamma104 1 Gamma126 1 Gamma128 1 Gamma800 1 Gamma151 1 Gamma152 1\n");
    fprintf (thisfile, "  Gamma130 1 Gamma132 1 Gamma44  1 Gamma53  1 Gamma49  1 Gamma804 1 Gamma805 1\n");
#elif defined USING_NBASE40
    fprintf (thisfile, "  Gamma70  1 Gamma77  1 Gamma78  1 Gamma802 1 Gamma803 1 Gamma93  1\n");
    fprintf (thisfile, "  Gamma94  1 Gamma103 1 Gamma104 1 Gamma126 1 Gamma128 1 Gamma800 1 Gamma151 1 Gamma152 1\n");
    fprintf (thisfile, "  Gamma130 1 Gamma132 1 Gamma44  1 Gamma53  1 Gamma49  1 Gamma804 1 Gamma805 1 Gamma801 1\n");
#endif
    if (!uconstrain) fprintf (thisfile, "  Gamma998 1\n");
  }
  if (p==0){
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_PSUM   1  0 ! print sum of strange decay nodes\n");
    fprintf (thisfile, "\n* Print Gamma(tau -> X-(S=1) nu)");
#if defined USING_NBASE31
    fprintf (thisfile, "\n*Gamma110 = Gamma10  + Gamma16   + Gamma23   + Gamma28  + Gamma35  + Gamma40 + Gamma85 + Gamma89 + Gamma128\n");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1     9  0");
#elif defined USING_NBASE39
    fprintf (thisfile, "\n*Gamma110 = Gamma10  + Gamma16   + Gamma23   + Gamma28  + Gamma35  + Gamma40 + Gamma802 + Gamma803 + Gamma128\n");
    fprintf (thisfile, "*         + Gamma151 + Gamma130  + Gamma132  + Gamma44  + Gamma53\n");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1     14 0");
#elif defined USING_NBASE40
    fprintf (thisfile, "\n*Gamma110 = Gamma10  + Gamma16   + Gamma23   + Gamma28  + Gamma35  + Gamma40 + Gamma802 + Gamma803 + Gamma128\n");
    fprintf (thisfile, "*         + Gamma151 + Gamma130  + Gamma132  + Gamma44  + Gamma53  + Gamma801\n");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1     15 0");
#endif
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_01  4   1 ! Gamma10");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_02  6   1 ! Gamma16");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_03  8   1 ! Gamma23");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_04  10  1 ! Gamma28");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_05  12  1 ! Gamma35");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_06  14  1 ! Gamma40");
#if defined USING_NBASE31
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_07  22  1 ! Gamma85");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_08  23  1 ! Gamma89");
#elif defined USING_NBASE39 || defined USING_NBASE40
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_07  22  1 ! Gamma802");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_08  23  1 ! Gamma803");
#endif
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_09  28  1 ! Gamma128");
#if defined USING_NBASE39 || defined USING_NBASE40
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_10  30  1 ! Gamma151");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_11  32  1 ! Gamma130");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_12  33  1 ! Gamma132");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_13  34  1 ! Gamma44");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_14  35  1 ! Gamma53");
#if defined USING_NBASE40
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_15  39  1 ! Gamma801");
#endif
#endif
    fprintf (thisfile, "\n");
  }
  if (p==1) {
    fprintf (thisfile, "\n* --- compute Gamma(tau -> Xs nu) / G(total)\n");
    fprintf (thisfile, "COMBOFQUANT Gamma110\n");
#if defined USING_NBASE31
    fprintf (thisfile, " 1 Gamma10  1 Gamma16  1 Gamma23  1 Gamma28  1 Gamma35  1 Gamma40  1 Gamma85  1 Gamma89  1 Gamma128\n");
#elif defined USING_NBASE39
    fprintf (thisfile, " 1 Gamma10  1 Gamma16  1 Gamma23  1 Gamma28  1 Gamma35  1 Gamma40  1 Gamma802 1 Gamma803 1 Gamma128\n");
    fprintf (thisfile, " 1 Gamma151 1 Gamma130 1 Gamma132 1 Gamma44  1 Gamma53\n");
#elif defined USING_NBASE40
    fprintf (thisfile, " 1 Gamma10  1 Gamma16  1 Gamma23  1 Gamma28  1 Gamma35  1 Gamma40  1 Gamma802 1 Gamma803 1 Gamma128\n");
    fprintf (thisfile, " 1 Gamma151 1 Gamma130 1 Gamma132 1 Gamma44  1 Gamma53  1 Gamma801\n");
#endif
  }
  fprintf (thisfile, "\nCALL CHI2_N_SYM\n");
  fprintf (thisfile, "\nEND\n");
}
// ----------------------------------------------------------------------
int main(int argc, char* argv[]){
  //Argument variables
  int uconstrain   = (argc>1) ? atoi(argv[1]) : 0; // 1: unitarity constrained; 0 : not constrained
  int doalephhcorr = (argc>2) ? atoi(argv[2]) : 1; // 1: do aleph hcorr; 0: dont
  //
  string sconstrain = (uconstrain == 1) ? "constrained" : "unconstrained";
  string salephhcorr = (doalephhcorr == 1) ? "_aleph_hcorr" : "";
  //
  gSystem->Load("libMatrix");
  //
  int ipar,iimeas,jjmeas,jpar;
  int i,j;
  int inode,jnode;
  int ibase,jbase;
  int im, jm;
  //
  const int weakcompare=1;
  //
  // READ BASE PARAMETERS
  //
#if defined USING_NBASE31
  const int nbase=31;
#elif defined USING_NBASE39
  const int nbase=39;
#elif defined USING_NBASE40
  const int nbase=40;
#endif
  int nbase_u = (uconstrain) ? nbase-1 : nbase;
  vector<int> basequan;
  vector<int> basegamma;
  vector<int> baseparm;
  double      baseseed[nbase];
  double      baseseed_orig[nbase];
  double      basefitval_orig[nbase];
  double      basefiterr_orig[nbase];
  double      basescalerr_orig[nbase];
  double      basescalfac_orig[nbase];
  string      basenodename[nbase];
  string      basetitle[nbase];
  //
  string basefilename = "base_def.txt";
  ifstream ifsbase(basefilename.data());
  if (!ifsbase.good()) {
    cout << Form("Cannot open input file : %s\n", basefilename.data()); 
    exit(1);
  } else {
    cout << Form("Read base definitions from : %s\n", basefilename.data());
  }
  //
  ibase=0;
  char buffer[256]; 
  while(ifsbase.good()) {
    if (ifsbase.eof()) break;
    char firstch(' ') ; ifsbase.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip such lines
      ifsbase.ignore(256,'\n') ;
    } else if (firstch=='*') {  // measurement line
      int dummy_int;       
      ifsbase >> dummy_int ; basequan.push_back(dummy_int);
      ifsbase >> dummy_int ; basegamma.push_back(dummy_int);
      ifsbase >> dummy_int ; baseparm.push_back(dummy_int);
      ifsbase >> basenodename[ibase];
      ifsbase >> baseseed_orig[ibase];
      ifsbase >> basefitval_orig[ibase];
      ifsbase >> basefiterr_orig[ibase];
      ifsbase >> basescalerr_orig[ibase];
      ifsbase >> basescalfac_orig[ibase];
      ifsbase.getline(buffer,256,'\n');
      string stringbuffer=string(buffer);
      basetitle[ibase]="";
      bool first=false;
      for (string::iterator si=stringbuffer.begin();si!=stringbuffer.end();++si){
	if (13 == (int)*si) break; // ^M is 13
	if (' '!=*si) first=true;
	if (first) {
	  if (' ' == *si && ' ' == *(si+1) ) break;
	  basetitle[ibase]+=*si;
	}
      }
      baseseed[ibase] = baseseed_orig[ibase];
      //      baseseed[ibase] = basefitval_orig[ibase];
      // cout << basequan[ibase] << " " << basegamma[ibase] << " " << baseparm[ibase] << " " << basenodename[ibase] << " " << baseseed[ibase] << " " << basetitle[ibase] << endl;
      ++ibase;
    }
  }
  //
  int          baseorder[nbase];
  string       baselatex[nbase];
  //
#if defined USING_NBASE31
  baseorder[M_GAMMA5  ] =  1; baselatex[M_GAMMA5  ] = "$e^- \\bar{\\nu}_e \\nu_\\tau$";
  baseorder[M_GAMMA3  ] =  2; baselatex[M_GAMMA3  ] = "$\\mu^- \\bar{\\nu}_\\mu \\nu_\\tau$";
  //
  baseorder[M_GAMMA9  ] =  3; baselatex[M_GAMMA9  ] = "$\\pi^- \\nu_\\tau$";
  baseorder[M_GAMMA14 ] =  4; baselatex[M_GAMMA14 ] = "$\\pi^- \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA20 ] =  5; baselatex[M_GAMMA20 ] = "$\\pi^- 2\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA27 ] =  6; baselatex[M_GAMMA27 ] = "$\\pi^- 3\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA30 ] =  7; baselatex[M_GAMMA30 ] = "$h^- 4\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\eta)$";
  baseorder[M_GAMMA37 ] =  8; baselatex[M_GAMMA37 ] = "$K^- K^0 \\nu_\\tau$";
  baseorder[M_GAMMA42 ] =  9; baselatex[M_GAMMA42 ] = "$K^- \\pi^0 K^0 \\nu_\\tau$";
  baseorder[M_GAMMA47 ] = 10; baselatex[M_GAMMA47 ] = "$\\pi^- K_S^0 K_S^0 \\nu_\\tau$";
  baseorder[M_GAMMA48 ] = 11; baselatex[M_GAMMA48 ] = "$\\pi^- K_S^0 K_L^0 \\nu_\\tau$";
  baseorder[M_GAMMA62 ] = 12; baselatex[M_GAMMA62 ] = "$\\pi^- \\pi^- \\pi^+ \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega)$";
  baseorder[M_GAMMA70 ] = 13; baselatex[M_GAMMA70 ] = "$\\pi^- \\pi^- \\pi^+ \\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega)$";
  baseorder[M_GAMMA77 ] = 14; baselatex[M_GAMMA77 ] = "$h^- h^- h^+ 2\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega,\\eta)$";
  baseorder[M_GAMMA78 ] = 15; baselatex[M_GAMMA78 ] = "$h^- h^- h^+ 3\\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA93 ] = 16; baselatex[M_GAMMA93 ] = "$\\pi^- K^- K^+ \\nu_\\tau$";
  baseorder[M_GAMMA94 ] = 17; baselatex[M_GAMMA94 ] = "$\\pi^- K^- K^+ \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA103] = 18; baselatex[M_GAMMA103] = "$3h^- 2h^+ \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA104] = 19; baselatex[M_GAMMA104] = "$3h^- 2h^+ \\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA126] = 20; baselatex[M_GAMMA126] = "$\\pi^- \\pi^0 \\eta \\nu_\\tau$";
  baseorder[M_GAMMA150] = 21; baselatex[M_GAMMA150] = "$h^- \\omega \\nu_\\tau$";
  baseorder[M_GAMMA152] = 22; baselatex[M_GAMMA152] = "$h^- \\pi^0 \\omega \\nu_\\tau$";
  //
  baseorder[M_GAMMA10 ] = 23; baselatex[M_GAMMA10 ] = "$K^- \\nu_\\tau$";
  baseorder[M_GAMMA16 ] = 24; baselatex[M_GAMMA16 ] = "$K^- \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA23 ] = 25; baselatex[M_GAMMA23 ] = "$K^- 2\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA28 ] = 26; baselatex[M_GAMMA28 ] = "$K^- 3\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\eta)$";
  baseorder[M_GAMMA35 ] = 27; baselatex[M_GAMMA35 ] = "$\\bar{K}^0 \\pi^- \\nu_\\tau$";
  baseorder[M_GAMMA40 ] = 28; baselatex[M_GAMMA40 ] = "$\\bar{K}^0 \\pi^- \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA85 ] = 29; baselatex[M_GAMMA85 ] = "$K^- \\pi^- \\pi^+ \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA89 ] = 30; baselatex[M_GAMMA89 ] = "$K^- \\pi^- \\pi^+ \\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\eta)$";
  baseorder[M_GAMMA128] = 31; baselatex[M_GAMMA128] = "$K^- \\eta \\nu_\\tau$";
  //
#elif defined USING_NBASE39 || defined USING_NBASE40
  //
  baseorder[M_GAMMA5  ] =  1; baselatex[M_GAMMA5  ] = "$e^- \\bar{\\nu}_e \\nu_\\tau$";
  baseorder[M_GAMMA3  ] =  2; baselatex[M_GAMMA3  ] = "$\\mu^- \\bar{\\nu}_\\mu \\nu_\\tau$";
  //
  baseorder[M_GAMMA9  ] =  3; baselatex[M_GAMMA9  ] = "$\\pi^- \\nu_\\tau$";
  baseorder[M_GAMMA14 ] =  4; baselatex[M_GAMMA14 ] = "$\\pi^- \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA20 ] =  5; baselatex[M_GAMMA20 ] = "$\\pi^- 2\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA27 ] =  6; baselatex[M_GAMMA27 ] = "$\\pi^- 3\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA30 ] =  7; baselatex[M_GAMMA30 ] = "$h^- 4\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\eta)$";
  baseorder[M_GAMMA37 ] =  8; baselatex[M_GAMMA37 ] = "$K^- K^0 \\nu_\\tau$";
  baseorder[M_GAMMA42 ] =  9; baselatex[M_GAMMA42 ] = "$K^- \\pi^0 K^0 \\nu_\\tau$";
  baseorder[M_GAMMA47 ] = 10; baselatex[M_GAMMA47 ] = "$\\pi^- K_S^0 K_S^0 \\nu_\\tau$";
  baseorder[M_GAMMA48 ] = 11; baselatex[M_GAMMA48 ] = "$\\pi^- K_S^0 K_L^0 \\nu_\\tau$";
  baseorder[M_GAMMA804] = 12; baselatex[M_GAMMA804] = "$\\pi^- K_L^0 K_L^0 \\nu_\\tau$";
  baseorder[M_GAMMA49 ] = 13; baselatex[M_GAMMA49 ] = "$\\pi^- \\pi^0 K^0 \\bar{K}^0 \\nu_\\tau$";
  baseorder[M_GAMMA62 ] = 14; baselatex[M_GAMMA62 ] = "$\\pi^- \\pi^- \\pi^+ \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega)$";
  baseorder[M_GAMMA70 ] = 15; baselatex[M_GAMMA70 ] = "$\\pi^- \\pi^- \\pi^+ \\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega)$";
  baseorder[M_GAMMA77 ] = 16; baselatex[M_GAMMA77 ] = "$h^- h^- h^+ 2\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega,\\eta)$";
  baseorder[M_GAMMA78 ] = 17; baselatex[M_GAMMA78 ] = "$h^- h^- h^+ 3\\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA93 ] = 18; baselatex[M_GAMMA93 ] = "$\\pi^- K^- K^+ \\nu_\\tau$";
  baseorder[M_GAMMA94 ] = 19; baselatex[M_GAMMA94 ] = "$\\pi^- K^- K^+ \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA103] = 20; baselatex[M_GAMMA103] = "$3h^- 2h^+ \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA104] = 21; baselatex[M_GAMMA104] = "$3h^- 2h^+ \\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA126] = 22; baselatex[M_GAMMA126] = "$\\pi^- \\pi^0 \\eta \\nu_\\tau$";
  baseorder[M_GAMMA800] = 23; baselatex[M_GAMMA800] = "$\\pi^- \\omega \\nu_\\tau$";
  baseorder[M_GAMMA152] = 24; baselatex[M_GAMMA152] = "$h^- \\pi^0 \\omega \\nu_\\tau$";
  baseorder[M_GAMMA805] = 25; baselatex[M_GAMMA805] = "$a_1^- (\\to \\pi^- \\gamma) \\nu_\\tau$";
  //
  baseorder[M_GAMMA10 ] = 26; baselatex[M_GAMMA10 ] = "$K^- \\nu_\\tau$";
  baseorder[M_GAMMA16 ] = 27; baselatex[M_GAMMA16 ] = "$K^- \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA23 ] = 28; baselatex[M_GAMMA23 ] = "$K^- 2\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0)$";
  baseorder[M_GAMMA28 ] = 29; baselatex[M_GAMMA28 ] = "$K^- 3\\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\eta)$";
  baseorder[M_GAMMA35 ] = 30; baselatex[M_GAMMA35 ] = "$\\bar{K}^0 \\pi^- \\nu_\\tau$";
  baseorder[M_GAMMA40 ] = 31; baselatex[M_GAMMA40 ] = "$\\bar{K}^0 \\pi^- \\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA44 ] = 32; baselatex[M_GAMMA44 ] = "$\\bar{K}^0 \\pi^- 2\\pi^0 \\nu_\\tau$";
  baseorder[M_GAMMA53 ] = 33; baselatex[M_GAMMA53 ] = "$\\bar{K}^0 h^- h^- h^+ \\nu_\\tau$";
  baseorder[M_GAMMA802] = 34; baselatex[M_GAMMA802] = "$K^- \\pi^- \\pi^+ \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega)$";
  baseorder[M_GAMMA803] = 35; baselatex[M_GAMMA803] = "$K^- \\pi^- \\pi^+ \\pi^0 \\nu_\\tau ~(\\mathrm{ex.~}K^0,\\omega,\\eta)$";
#if defined USING_NBASE39
  baseorder[M_GAMMA128] = 36; baselatex[M_GAMMA128] = "$K^- \\eta \\nu_\\tau$";
  baseorder[M_GAMMA130] = 37; baselatex[M_GAMMA130] = "$K^- \\pi^0 \\eta \\nu_\\tau$";
  baseorder[M_GAMMA132] = 38; baselatex[M_GAMMA132] = "$\\bar{K}^0 \\pi^- \\eta \\nu_\\tau$";
  baseorder[M_GAMMA151] = 39; baselatex[M_GAMMA151] = "$K^- \\omega \\nu_\\tau$";
#elif defined USING_NBASE40
  baseorder[M_GAMMA801] = 36; baselatex[M_GAMMA801] = "$K^- \\phi \\nu_\\tau (\\phi \\to KK)$";
  baseorder[M_GAMMA128] = 37; baselatex[M_GAMMA128] = "$K^- \\eta \\nu_\\tau$";
  baseorder[M_GAMMA130] = 38; baselatex[M_GAMMA130] = "$K^- \\pi^0 \\eta \\nu_\\tau$";
  baseorder[M_GAMMA132] = 39; baselatex[M_GAMMA132] = "$\\bar{K}^0 \\pi^- \\eta \\nu_\\tau$";
  baseorder[M_GAMMA151] = 40; baselatex[M_GAMMA151] = "$K^- \\omega \\nu_\\tau$";
#endif
#endif
  //
  for (i=0;i<nbase;++i) {
    for (ibase=0;ibase<nbase;++ibase) {
      if (baseorder[ibase]==i) {
	cout << "ibase = " << ibase << endl;
	cout << "title = " << basetitle[ibase] << endl;
	cout << "latex = " << baselatex[ibase] << endl << endl;
      }
    }
  }
  //
  // DEFINE COEFFICIENTS
  //
  const int ncoef=32;
  string coefnames[ncoef];
  coefnames[C_ETA3PIZUSAGE1    ] = "B(eta -> pi0 pi0 pi0)";
  coefnames[C_ETA3PIZUSAGE2    ] = "B(eta -> pi0 pi0 pi0)"; // [used for some older usage of pi- pi0 eta nu(tau)]
  coefnames[C_ETA3PIZUSAGE3    ] = "B(eta -> pi0 pi0 pi0)";//  [used for some older usage of K- eta nu(tau) and pi- pi0 eta nu(tau)]
  coefnames[C_ETAPIMPIPPIZ     ] = "B(eta -> pi- pi+ pi0)";
  coefnames[C_ETA3PIALL        ] = "B(eta -> 3pi0, pi-pi+pi0)"; // B(eta -> pi0 pi0 pi0) + B(eta -> pi- pi+ pi0) 
  coefnames[C_ETANEUTRALMODES1 ] = "B(eta -> neutrals)"; // sum of neutral modes, mostly = B(eta -> pi0 pi0 pi0) + B(eta -> gamma gamma) [used for K- eta nu(tau)]
  coefnames[C_ETANEUTRALMODES2 ] = "B(eta -> neutrals)"; // sum of neutral modes, mostly = B(eta -> pi0 pi0 pi0) + B(eta -> gamma gamma) [used for pi- pi0 eta nu(tau)]
  coefnames[C_ETACHARGEDMODES  ] = "B(eta -> charged)"; // sum of charged modes, mostly = B(eta -> pi- pi+ pi0) + B(eta -> pi- pi+ gamma)
  coefnames[C_KS2PIZ           ] = "B(KS -> pi0 pi0)"; // 
  coefnames[C_KS2PIZby2        ] = "0.5 * B(KS -> pi0 pi0)";
  coefnames[C_KS2PIZby2PlusHalf] = "(0.5 * B(KS -> pi0 pi0) + 0.5)";
  coefnames[C_KS2PIZSQUARED    ] = "B(KS -> pi0 pi0) * B(KS -> pi0 pi0)";
  coefnames[C_KS2PIZSQUAREPlus1] = "(B(KS -> pi0 pi0) * B(KS -> pi0 pi0) + 1)";
  coefnames[C_KS2PI            ] = "B(KS -> pi- pi+)";
  coefnames[C_KS2PIby2         ] = "0.5 * B(KS -> pi- pi+)";
  coefnames[C_KS2PI_X_2PIZ_X_2 ] = "2 * B(KS -> pi- pi+) * B(KS -> pi- pi0)";
  coefnames[C_OMPIMPIPPIZ      ] = "B(omega -> pi- pi+ pi0)";
  coefnames[C_OMPIMPIP         ] = "B(omega -> pi- pi+)";
  coefnames[C_OM3PIPlus2PI     ] = "B(omega -> pi- pi+ pi0) + B(omega -> pi- pi+)";
  coefnames[C_OM2PIZGAMMA      ] = "B(omega -> pi0 gamma)";
  coefnames[C_PHI2KMKP         ] = "B(phi->K-K+)/(B(phi->K-K+)+B(phi->KSKL))";
  coefnames[C_PHI2KSKL         ] = "B(phi->KSKL)/(B(phi->K-K+)+B(phi->KSKL))"; 
  coefnames[C_PHI2KK_KS2PI     ] = "(B(phi->K-K+)+B(phi->KSKL)*B(KS->pi-pi+))/(B(phi->K-K+)+B(phi->KSKL))"; 
  coefnames[C_PHI2KSKL_KS2PI   ] = "(B(phi->KSKL)*B(KS->pi-pi+))/(B(phi->K-K+)+B(phi->KSKL))"; 
  coefnames[C_PHI2KSKL_KS2PIZ  ] = "(B(phi->KSKL)*B(KS->pi0pi0))/(B(phi->K-K+)+B(phi->KSKL))"; 
  coefnames[C_Gam132_3PI3PIZ_KL] = "(0.5*(B(eta->pi-pi+pi0) + B(KS->2pi0)*B(eta->pi-pi+pi0) + B(KS->pi-pi+)*B(eta->3pi0)))";
  coefnames[C_Gam132_ETANeuby2 ] = "0.5 * B(eta -> neutrals)";
  coefnames[C_ONE              ] = "ONE";
  coefnames[C_HALF             ] = "HALF";
  coefnames[C_TWO              ] = "TWO";
  coefnames[C_ONETHIRD         ] = "ONETHIRD";
  coefnames[C_TWOTHIRD         ] = "TWOTHIRD";
  //
  double coef[ncoef];
#if defined USING_NBASE31
  // original values used in pdgfit/s035-fit-no-babar-belle.fit
  coef[C_ETA3PIZUSAGE1    ] = 3.2200E-01; // B(eta -> pi0 pi0 pi0)
  coef[C_ETA3PIZUSAGE2    ] = 3.1900E-01; // B(eta -> pi0 pi0 pi0) [used for some older usage of pi- pi0 eta nu(tau)]
  coef[C_ETA3PIZUSAGE3    ] = 3.2500E-01; // B(eta -> pi0 pi0 pi0) [used for some older usage of K- eta nu(tau) and pi- pi0 eta nu(tau)]
  coef[C_ETAPIMPIPPIZ     ] = 2.2600E-01; // B(eta -> pi- pi+ pi0)
  coef[C_ETA3PIALL        ] = 5.5300E-01; // B(eta -> pi0 pi0 pi0) + B(eta -> pi- pi+ pi0)
  coef[C_ETANEUTRALMODES1 ] = 7.1500E-01; // sum of neutral modes, mostly = B(eta -> pi0 pi0 pi0) + B(eta -> gamma gamma) [used for K- eta nu(tau)]
  coef[C_ETANEUTRALMODES2 ] = 7.0800E-01; // sum of neutral modes, mostly = B(eta -> pi0 pi0 pi0) + B(eta -> gamma gamma) [used for pi- pi0 eta nu(tau)]
  coef[C_ETACHARGEDMODES  ] = 2.8500E-01; // sum of charged modes, mostly = B(eta -> pi- pi+ pi0) + B(eta -> pi- pi+ gamma)
  coef[C_KS2PIZ           ] = 3.1390E-01; // B(KS -> pi0 pi0)
  coef[C_KS2PIZby2        ] = 1.5700E-01; // 0.5 x B(KS -> pi0 pi0)
  coef[C_KS2PIZby2PlusHalf] = 6.5690E-01; // 0.5 x B(KS -> pi0 pi0) + 0.5
  coef[C_KS2PIZSQUARED    ] = 9.8500E-02; // B(KS -> pi0 pi0) * B(KS -> pi0 pi0)
  coef[C_KS2PIZSQUAREPlus1] = 1.0985E+00; // B(KS -> pi0 pi0) * B(KS -> pi0 pi0) + 1
  coef[C_KS2PIby2         ] = 3.4310E-01; // 0.5 x B(KS -> pi- pi+)
  coef[C_KS2PI            ] = 6.8610E-01; //       B(KS -> pi- pi+)
  coef[C_KS2PI_X_2PIZ_X_2 ] = 4.3070E-01; //   2 x B(KS -> pi- pi+) + B(KS -> pi0 pi0)
  coef[C_OMPIMPIPPIZ      ] = 8.8800E-01; // B(omega -> pi- pi+ pi0)
  coef[C_OMPIMPIP         ] = 1.7000E-02; // B(omega -> pi- pi+)
  coef[C_OM3PIPlus2PI     ] = 9.1010E-01; // B(omega -> pi- pi+ pi0) + B(omega -> pi- pi+)
  coef[C_OM2PIZGAMMA      ] = 9.0000E-02; // B(omega -> pi0 gamma)
  coef[C_PHI2KMKP         ] = 0.588448;   // B(phi->K-K+)/(B(phi->K-K+)+B(phi->KSKL)) using B(K-K+)=0.489+-0.005, B(KSKL)=0.342+-0.004 [PDG09]
  coef[C_PHI2KSKL         ] = 0.411552;   // B(phi->KSKL)/(B(phi->K-K+)+B(phi->KSKL)) using B(K-K+)=0.489+-0.005, B(KSKL)=0.342+-0.004 [PDG09]
  coef[C_PHI2KK_KS2PI     ] = 0.870814;   // ( .489 + .342*6.8610E-01 ) / ( .489+.342 ) = 0.870814
  coef[C_PHI2KSKL_KS2PI   ] = 0.282366;   // ( .342*6.8610E-01 ) / ( .489+.342 ) = 0.282366
  coef[C_PHI2KSKL_KS2PIZ  ] = 0.129186;   // ( .342*3.1390E-01 ) / ( .489+.342 ) = 0.129186
  coef[C_Gam132_3PI3PIZ_KL] = 0.261287;   // 0.5 * ( 22.74e-2 + 30.69e-2 * 22.74e-2 + 69.20e-2 * 32.57e-2) = 0.261287 //Kbar0 pi- eta nu(tau) -> 3pi- 3pi0, 3pi- pi0 KL0
  coef[C_Gam132_ETANeuby2 ] = 0.3595;     // 0.5 * BR_eta_neutral
  coef[C_ONE              ] = 1.0000E+00;
  coef[C_HALF             ] = 5.0000E-01;
  coef[C_TWO              ] = 2.0000E+00;
  coef[C_ONETHIRD         ] = 1.0E+00/3.0E+00;
  coef[C_TWOTHIRD         ] = 2.0E+00/3.0E+00;
#elif defined USING_NBASE39 || defined USING_NBASE40
  // PDG2010
  double BR_eta_2gamgam     = 39.31e-2; // ( 39.31 +- 0.20 ) %
  double BR_eta_neutral     = 71.90e-2; // ( 71.90 +- 0.34 ) %
  double BR_eta_3piz        = 32.57e-2; // ( 32.57 +- 0.23 ) %
  double BR_eta_pimpippiz   = 22.74e-2; // ( 22.74 +- 0.28 ) %
  double BR_eta_charged     = 28.10e-2; // ( 28.10 +- 0.34 ) %
  double BR_KS_2piz         = 30.69e-2; // ( 30.69 +- 0.05 ) %
  double BR_KS_pimpip       = 69.20e-2; // ( 69.20 +- 0.05 ) %
  double BR_om_pimpippiz    = 89.2e-2;  // ( 89.2  +- 0.7 ) %
  double BR_om_pimpip       = 1.53e-2;  // ( 1.53  +  0.11 - 0.13 ) %
  double BR_om_pizgamma     = 8.28e-2;  // ( 8.28  +- 0.28 ) %
  double BR_phi_KmKp        = 48.9e-2;  // ( 48.9  +- 0.5 ) %
  double BR_phi_KSKL        = 34.2e-2;  // ( 34.2  +- 0.4 ) %
  //
  coef[C_ETA3PIZUSAGE1    ] = BR_eta_3piz; // B(eta -> pi0 pi0 pi0)
  coef[C_ETA3PIZUSAGE2    ] = BR_eta_3piz; // B(eta -> pi0 pi0 pi0) [used for some older usage of pi- pi0 eta nu(tau)]
  coef[C_ETA3PIZUSAGE3    ] = BR_eta_3piz; // B(eta -> pi0 pi0 pi0) [used for some older usage of K- eta nu(tau) and pi- pi0 eta nu(tau)]
  coef[C_ETAPIMPIPPIZ     ] = BR_eta_pimpippiz; // B(eta -> pi- pi+ pi0)
  coef[C_ETA3PIALL        ] = BR_eta_3piz + BR_eta_pimpippiz; //B(eta -> pi0 pi0 pi0) + B(eta -> pi- pi+ pi0)
  coef[C_ETANEUTRALMODES1 ] = BR_eta_neutral; // sum of neutral modes, mostly = B(eta -> pi0 pi0 pi0) + B(eta -> gamma gamma) [used for K- eta nu(tau)]
  coef[C_ETANEUTRALMODES2 ] = BR_eta_neutral; // sum of neutral modes, mostly = B(eta -> pi0 pi0 pi0) + B(eta -> gamma gamma) [used for pi- pi0 eta nu(tau)]
  coef[C_ETACHARGEDMODES  ] = BR_eta_charged; // sum of charged modes, mostly = B(eta -> pi- pi+ pi0) + B(eta -> pi- pi+ gamma)
  coef[C_KS2PIZ           ] = BR_KS_2piz;  // B(KS -> pi0 pi0)
  coef[C_KS2PIZby2        ] = 0.5*BR_KS_2piz; // 0.5 x B(KS -> pi0 pi0)
  coef[C_KS2PIZby2PlusHalf] = (0.5*BR_KS_2piz) + 0.5; // 0.5 x B(KS -> pi0 pi0) + 0.5
  coef[C_KS2PIZSQUARED    ] = BR_KS_2piz*BR_KS_2piz; // B(KS -> pi0 pi0) * B(KS -> pi0 pi0)
  coef[C_KS2PIZSQUAREPlus1] = (BR_KS_2piz*BR_KS_2piz)+1; // B(KS -> pi0 pi0) * B(KS -> pi0 pi0) + 1
  coef[C_KS2PI            ] = BR_KS_pimpip; // B(KS -> pi- pi+)
  coef[C_KS2PIby2         ] = 0.5*BR_KS_pimpip; // 0.5 x B(KS -> pi- pi+)
  coef[C_KS2PI_X_2PIZ_X_2 ] = 2*BR_KS_pimpip*BR_KS_2piz; // 2 * B(KS -> pi- pi+) + B(KS -> pi0 pi0)
  coef[C_OMPIMPIPPIZ      ] = BR_om_pimpippiz; // B(omega -> pi- pi+ pi0)
  coef[C_OMPIMPIP         ] = BR_om_pimpip; // B(omega -> pi- pi+)
  coef[C_OM3PIPlus2PI     ] = BR_om_pimpippiz + BR_om_pimpip; // B(omega -> pi- pi+ pi0) + B(omega -> pi- pi+)
  coef[C_OM2PIZGAMMA      ] = BR_om_pizgamma; // B(omega -> pi0 gamma)
  coef[C_PHI2KMKP         ] = BR_phi_KmKp/(BR_phi_KmKp+BR_phi_KSKL);// B(phi->K-K+)/(B(phi->K-K+)+B(phi->KSKL))
  coef[C_PHI2KSKL         ] = BR_phi_KSKL/(BR_phi_KmKp+BR_phi_KSKL);// B(phi->KSKL)/(B(phi->K-K+)+B(phi->KSKL))
  coef[C_PHI2KK_KS2PI     ] = (BR_phi_KmKp + BR_phi_KSKL*BR_KS_pimpip)/(BR_phi_KmKp+BR_phi_KSKL);// (B(phi->K-K+) + B(phi->KSKL)*B(KS->pi-pi+))/(B(phi->K-K+)+B(phi->KSKL))
  coef[C_PHI2KSKL_KS2PI   ] = (BR_phi_KSKL*BR_KS_pimpip)/(BR_phi_KmKp+BR_phi_KSKL);// (B(phi->KSKL)*B(KS->pi-pi+))/(B(phi->K-K+)+B(phi->KSKL))
  coef[C_PHI2KSKL_KS2PIZ  ] = (BR_phi_KSKL*BR_KS_2piz)/(BR_phi_KmKp+BR_phi_KSKL);// (B(phi->KSKL)*B(KS->pi0pi0))/(B(phi->K-K+)+B(phi->KSKL))
  coef[C_Gam132_3PI3PIZ_KL] = 0.5 * BR_eta_pimpippiz + 0.5 * BR_KS_2piz * BR_eta_pimpippiz + 0.5 * BR_KS_pimpip * BR_eta_3piz;//Kbar0 pi- eta nu(tau) -> 3pi- 3pi0, 3pi- pi0 KL0
  coef[C_Gam132_ETANeuby2 ] = 0.5 * BR_eta_neutral;
  coef[C_ONE              ] = 1.0000E+00;
  coef[C_HALF             ] = 5.0000E-01;
  coef[C_TWO              ] = 2.0000E+00;
  coef[C_ONETHIRD         ] = 1.0E+00/3.0E+00;
  coef[C_TWOTHIRD         ] = 2.0E+00/3.0E+00;
#endif
  //
  // READ INPUT NODES
  // 
#if defined USING_NBASE31
  const int nnode=99;
#elif defined USING_NBASE39 
  const int nnode=113; 
#elif defined USING_NBASE40
  const int nnode=115; 
#endif
  vector<string> nodegammaname;
  vector<string> nodename;
  vector<int>    node_num_parm[nnode];
  vector<double> node_num_coef[nnode];
  vector<int>    node_den_parm[nnode];
  vector<double> node_den_coef[nnode];
  string         nodetitle[nnode];
  //
  inode=0;//0
  nodegammaname.push_back("Gamma128");     
  nodename.push_back("S035B20");    
  nodetitle[inode]="G(K- eta nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- eta nu(tau) :: Gamma128
  ++inode;//1
  nodegammaname.push_back("Gamma19by13");  
  nodename.push_back("S035B21"); 
  nodetitle[inode]="G(h- 2pi0 nu(tau) (ex. K0)) / G(h- pi0 nu(tau))";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_den_parm[inode].push_back(16 );       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_den_parm[inode].push_back(182);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  ++inode;//2
  nodegammaname.push_back("Gamma26by13"); 
  nodename.push_back("S035B22"); 
  nodetitle[inode]="G(h- 3pi0 nu(tau)) / G(h- pi0 nu(tau))";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE1    ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_den_parm[inode].push_back(16 );       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_den_parm[inode].push_back(182);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  ++inode;//3
  nodegammaname.push_back("Gamma30"); 
  nodename.push_back("S035B23"); 
  nodetitle[inode]="G(h- 4pi0 nu(tau) (ex. K0, eta)) / G(total)";
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  ++inode;//4
  nodegammaname.push_back("Gamma76by54"); 
  nodename.push_back("S035B25"); 
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex. K0)) / G(h- h- h+ >=0 neutrals >=0 K(L)0 nu(tau))";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(117);       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_den_parm[inode].push_back(118);       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_den_parm[inode].push_back(119);       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_den_parm[inode].push_back(204);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_den_parm[inode].push_back(214);       node_den_coef[inode].push_back(coef[C_KS2PI_X_2PIZ_X_2 ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_den_parm[inode].push_back(240);       node_den_coef[inode].push_back(coef[C_KS2PI            ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_den_parm[inode].push_back(5  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_den_parm[inode].push_back(62 );       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE31
  node_den_parm[inode].push_back(260);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(802);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_den_parm[inode].push_back(803);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// pi- omega nu(tau) :: Gamma800
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
  node_den_parm[inode].push_back(266);       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- pi0 eta nu(tau) :: Gamma130
  node_den_parm[inode].push_back(267);       node_den_coef[inode].push_back(coef[C_Gam132_3PI3PIZ_KL]);// Kbar0 pi- eta nu(tau) -> 3pi- 3pi0, 3pi- pi0 KL0 :: Gamma132 
  node_den_parm[inode].push_back(244);       node_den_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
#if defined USING_NBASE40
  node_den_parm[inode].push_back(801);       node_den_coef[inode].push_back(coef[C_PHI2KK_KS2PI     ]);// K- phi nu(tau) (phi->K-K+, KS(->pi-pi+)KL) :: Gamma801
#endif
#endif
  ++inode;//5
  nodegammaname.push_back("Gamma152by54"); 
  nodename.push_back("S035B26"); 
  nodetitle[inode]="G(h- omega pi0 nu(tau)) / G(h- h- h+ >=0 neutrals >=0 K(L)0 nu(tau))";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(117);       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_den_parm[inode].push_back(118);       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_den_parm[inode].push_back(119);       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_den_parm[inode].push_back(204);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_den_parm[inode].push_back(214);       node_den_coef[inode].push_back(coef[C_KS2PI_X_2PIZ_X_2 ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_den_parm[inode].push_back(240);       node_den_coef[inode].push_back(coef[C_KS2PI            ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_den_parm[inode].push_back(5  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_den_parm[inode].push_back(62 );       node_den_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE31
  node_den_parm[inode].push_back(260);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(802);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_den_parm[inode].push_back(803);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// pi- omega nu(tau) :: Gamma800
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
  node_den_parm[inode].push_back(266);       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- pi0 eta nu(tau) :: Gamma130
  node_den_parm[inode].push_back(267);       node_den_coef[inode].push_back(coef[C_Gam132_3PI3PIZ_KL]);// Kbar0 pi- eta nu(tau) -> 3pi- 3pi0, 3pi- pi0 KL0 :: Gamma132 
  node_den_parm[inode].push_back(244);       node_den_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
#if defined USING_NBASE40
  node_den_parm[inode].push_back(801);       node_den_coef[inode].push_back(coef[C_PHI2KK_KS2PI     ]);// K- phi nu(tau) (phi->K-K+, KS(->pi-pi+)KL) :: Gamma801
#endif
#endif
  ++inode;//6
  nodegammaname.push_back("Gamma152by76"); 
  nodename.push_back("S035B27"); 
  nodetitle[inode]="G(h- omega pi0 nu(tau)) / G(h- h- h+ 2pi0 nu(tau) (ex. K0))";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(266);       node_den_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//7
  nodegammaname.push_back("Gamma16"); 
  nodename.push_back("S035B29"); 
  nodetitle[inode]="G(K- pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  ++inode;//8
  nodegammaname.push_back("Gamma23"); 
  nodename.push_back("S035B30"); 
  nodetitle[inode]="G(K- 2pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  ++inode;//9
  nodegammaname.push_back("Gamma28"); 
  nodename.push_back("S035B31"); 
  nodetitle[inode]="G(K- 3pi0 nu(tau) (ex. K0, eta)) / G(total)";
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  ++inode;//10
  nodegammaname.push_back("Gamma35"); 
  nodename.push_back("S035B32"); 
  nodetitle[inode]="G(Kbar0 pi- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- nu(tau) :: Gamma35
  ++inode;//11
  nodegammaname.push_back("Gamma40"); 
  nodename.push_back("S035B33"); 
  nodetitle[inode]="G(Kbar0 pi- pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  ++inode;//12
  nodegammaname.push_back("Gamma42"); 
  nodename.push_back("S035B34"); 
  nodetitle[inode]="G(K- pi0 K0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  ++inode;//13
  nodegammaname.push_back("Gamma92"); 
  nodename.push_back("S035B37"); 
  nodetitle[inode]="G(pi- K- K+ >=0 neutrals nu(tau)) / G(total)";
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  ++inode;//14
  nodegammaname.push_back("Gamma33"); 
  nodename.push_back("S035B43"); 
  nodetitle[inode]="G(K(S)0 (particles)- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_HALF             ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_HALF             ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_HALF             ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_HALF             ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_Gam132_ETANeuby2 ]);// Kbar0 pi- eta nu(tau) :: Gamma132
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_HALF             ]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KSKL         ]);// K- phi nu(tau) (phi->K(S)0 K(L)0) :: Gamma801
#endif
#endif
  ++inode;//15
  nodegammaname.push_back("Gamma106"); 
  nodename.push_back("S035B45"); 
  nodetitle[inode]="G((5pi)- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ nu(tau) (ex. K0) :: Gamma103
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETA3PIALL        ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETA3PIALL        ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
#endif
  ++inode;//16
  nodegammaname.push_back("Gamma46"); 
  nodename.push_back("S035B51"); 
  nodetitle[inode]="G(pi- K0 Kbar0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
#if defined USING_NBASE31
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_TWO              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(804);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(L)0 K(L)0 nu(tau) :: Gamma804  
#endif
  ++inode;//17
  nodegammaname.push_back("Gamma66"); 
  nodename.push_back("S035B53"); 
  nodetitle[inode]="G(h- h- h+ pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//18
  nodegammaname.push_back("Gamma67"); 
  nodename.push_back("S035B54"); 
  nodetitle[inode]="G(h- h- h+ pi0 nu(tau) (ex. K0, omega)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//19
  nodegammaname.push_back("Gamma20"); 
  nodename.push_back("S035B55"); 
  nodetitle[inode]="G(pi- 2pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  ++inode;//20
  nodegammaname.push_back("Gamma27"); 
  nodename.push_back("S035B56"); 
  nodetitle[inode]="G(pi- 3pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  ++inode;//21
  nodegammaname.push_back("Gamma78"); 
  nodename.push_back("S035B57"); 
  nodetitle[inode]="G(h- h- h+ 3pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  ++inode;//22
  nodegammaname.push_back("Gamma152"); 
  nodename.push_back("S035B58"); 
  nodetitle[inode]="G(h- pi0 omega nu(tau)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- pi0 omega nu(tau) :: Gamma152
  ++inode;//23
  nodegammaname.push_back("Gamma76"); 
  nodename.push_back("S035B59"); 
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//24
  nodegammaname.push_back("Gamma57"); 
  nodename.push_back("S035B62"); 
  nodetitle[inode]="G(h- h- h+ nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
#endif
#endif
  ++inode;//25
  nodegammaname.push_back("Gamma55"); 
  nodename.push_back("S035B63"); 
  nodetitle[inode]="G(h- h- h+ >=0 neutrals nu(tau) (ex. K(S)0 --> pi- pi+) (``3-prong'')) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- pi0 eta nu(tau) :: Gamma130
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
#endif
#endif
  ++inode;//26
  nodegammaname.push_back("Gamma57by55"); 
  nodename.push_back("S035B64"); 
  nodetitle[inode]="G(h- h- h+ nu(tau) (ex. K0)) / G(h- h- h+ >=0 neutrals nu(tau) (ex. K(S)0 --> pi- pi+) (``3-prong''))";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
#endif
#endif
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(204);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_den_parm[inode].push_back(5  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE31
  node_den_parm[inode].push_back(260);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(802);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_den_parm[inode].push_back(803);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// pi- omega nu(tau) :: Gamma800
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
  node_den_parm[inode].push_back(266);       node_den_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- pi0 eta nu(tau) :: Gamma130
#if defined USING_NBASE40
  node_den_parm[inode].push_back(801);       node_den_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
#endif
#endif
  ++inode;//27
  nodegammaname.push_back("Gamma34"); 
  nodename.push_back("S035B67"); 
  nodetitle[inode]="G(Kbar0 h- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//28
  nodegammaname.push_back("Gamma39"); 
  nodename.push_back("S035B68"); 
  nodetitle[inode]="G(Kbar0 h- pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  ++inode;//29
  nodegammaname.push_back("Gamma47"); 
  nodename.push_back("S035B69"); 
  nodetitle[inode]="G(pi- K(S)0 K(S)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  ++inode;//30
  nodegammaname.push_back("Gamma58"); 
  nodename.push_back("S035B71"); 
  nodetitle[inode]="G(h- h- h+ nu(tau) (ex. K0, omega)) / G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
#endif
#endif
  ++inode;//31
  nodegammaname.push_back("Gamma77");
  nodename.push_back("S035B72"); 
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta)) / G(total)";
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  ++inode;//32
  nodegammaname.push_back("Gamma8"); 
  nodename.push_back("S035B73"); 
  nodetitle[inode]="G(h- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  ++inode;//33
  nodegammaname.push_back("Gamma18");
  nodename.push_back("S035B74"); 
  nodetitle[inode]="G(h- 2pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//34
  nodegammaname.push_back("Gamma1");
  nodename.push_back("S035B75"); 
  nodetitle[inode]="G((particles)- >=0 neutrals >=0 K0 nu(tau) (``1-prong'')) / G(total)";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// mu- nubar(mu) nu(tau) :: Gamma3
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES1 ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
#if defined USING_NBASE31
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_TWO              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(804);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(L)0 K(L)0 nu(tau) :: Gamma804  
#endif
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES2 ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES2 ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_Gam132_ETANeuby2 ]);// Kbar0 pi- eta nu(tau) :: Gamma132 
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_HALF             ]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KSKL_KS2PIZ  ]);// K- phi nu(tau) (phi->KS(->pi0pi0)KL) :: Gamma801
#endif
#endif
  ++inode;//35
  nodegammaname.push_back("Gamma65");
  nodename.push_back("S035B76"); 
  nodetitle[inode]="G(h- h- h+ pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//36
  nodegammaname.push_back("Gamma75");
  nodename.push_back("S035B77");
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PI_X_2PIZ_X_2 ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//37
  nodegammaname.push_back("Gamma64");
  nodename.push_back("S035B78"); 
  nodetitle[inode]="G(h- h- h+ >=1 pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//38
  nodegammaname.push_back("Gamma29"); 
  nodename.push_back("S035B79"); 
  nodetitle[inode]="G(h- 4pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE2    ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE2    ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//39
  nodegammaname.push_back("Gamma8by5");
  nodename.push_back("S035B97"); 
  nodetitle[inode]="G(h- nu(tau)) / G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  ++inode;//40
  nodegammaname.push_back("Gamma12");
  nodename.push_back("S035C01"); 
  nodetitle[inode]="G(h- >= 1pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE3    ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE3    ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE3    ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//41
  nodegammaname.push_back("Gamma25"); 
  nodename.push_back("S035C02"); 
  nodetitle[inode]="G(h- >= 3pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE3    ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE3    ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE3    ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//42
  nodegammaname.push_back("Gamma74");
  nodename.push_back("S035C03");
  nodetitle[inode]="G(h- h- h+ >= 2pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//43
  nodegammaname.push_back("Gamma48");
  nodename.push_back("S035C1");
  nodetitle[inode]="G(pi- K(S)0 K(L)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  ++inode;//44
  nodegammaname.push_back("Gamma59");
  nodename.push_back("S035C18"); 
  nodetitle[inode]="G(pi- pi- pi+ nu(tau)) / G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//45
  nodegammaname.push_back("Gamma60");
  nodename.push_back("S035C19");
  nodetitle[inode]="G(pi- pi- pi+ nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//46
  nodegammaname.push_back("Gamma62");
  nodename.push_back("S035C20");
  nodetitle[inode]="G(pi- pi- pi+ nu(tau) (ex. K0, omega)) / G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  ++inode;//47
  nodegammaname.push_back("Gamma85");
  nodename.push_back("S035C21");
  nodetitle[inode]="G(K- pi- pi+ nu(tau) (ex. K0)) / G(total)";
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//48
  nodegammaname.push_back("Gamma68");
  nodename.push_back("S035C22");
  nodetitle[inode]="G(pi- pi- pi+ pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//49
  nodegammaname.push_back("Gamma69");
  nodename.push_back("S035C23");
  nodetitle[inode]="G(pi- pi- pi+ pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//50
  nodegammaname.push_back("Gamma70");
  nodename.push_back("S035C24");
  nodetitle[inode]="G(pi- pi- pi+ pi0 nu(tau) (ex. K0, omega)) / G(total)";
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  ++inode;//51
  nodegammaname.push_back("Gamma88"); 
  nodename.push_back("S035C25");
  nodetitle[inode]="G(K- pi- pi+ pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//52
  nodegammaname.push_back("Gamma80"); 
  nodename.push_back("S035C31"); 
  nodetitle[inode]="G(K- pi- h+ nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//53
  nodegammaname.push_back("Gamma80by60"); 
  nodename.push_back("S035C32"); 
  nodetitle[inode]="G(K- pi- h+ nu(tau) (ex. K0)) / G(pi- pi- pi+ nu(tau) (ex. K0))";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#endif
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
#if defined USING_NBASE31
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//54
  nodegammaname.push_back("Gamma81"); 
  nodename.push_back("S035C33"); 
  nodetitle[inode]="G(K- pi- h+ pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//55
  nodegammaname.push_back("Gamma81by69");
  nodename.push_back("S035C34");
  nodetitle[inode]="G(K- pi- h+ pi0 nu(tau) (ex. K0)) / G(pi- pi- pi+ pi0 nu(tau) (ex. K0))";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//56
  nodegammaname.push_back("Gamma93by60");
  nodename.push_back("S035C35");
  nodetitle[inode]="G(pi- K- K+ nu(tau)) / G(pi- pi- pi+ nu(tau) (ex. K0))";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
#if defined USING_NBASE31
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//57
  nodegammaname.push_back("Gamma94by69");
  nodename.push_back("S035C36");
  nodetitle[inode]="G(pi- K- K+ pi0 nu(tau)) / G(pi- pi- pi+ pi0 nu(tau) (ex. K0))";
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
#endif
  ++inode;//58
  nodegammaname.push_back("Gamma38");
  nodename.push_back("S035C38");
  nodetitle[inode]="G(K- K0 >=0 pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//59
  nodegammaname.push_back("Gamma83");
  nodename.push_back("S035C40"); 
  nodetitle[inode]="G(K- pi- pi+ >=0 pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//60
  nodegammaname.push_back("Gamma110");
  nodename.push_back("S035C47"); 
  nodetitle[inode]="G(strange) / G(total)";
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- eta nu(tau) :: Gamma128
#if defined USING_NBASE31 
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- eta nu(tau) :: Gamma132
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- phi nu(tau) (phi->KK) :: Gamma801
#endif
#endif
  ++inode;//61
  nodegammaname.push_back("Gamma89");
  nodename.push_back("S035C54");
  nodetitle[inode]="G(K- pi- pi+ pi0 nu(tau) (ex. K0, eta)) / G(total)";
#if defined USING_NBASE31 
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//62
  nodegammaname.push_back("Gamma84");
  nodename.push_back("S035C6"); 
  nodetitle[inode]="G(K- pi- pi+ nu(tau)) / G(total)";
#if defined USING_NBASE31 
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#endif
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//63
  nodegammaname.push_back("Gamma87");
  nodename.push_back("S035C7");
  nodetitle[inode]="G(K- pi- pi+ pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
#if defined USING_NBASE31 
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//64
  nodegammaname.push_back("Gamma94"); 
  nodename.push_back("S035C8");
  nodetitle[inode]="G(pi- K- K+ pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  ++inode;//65
  nodegammaname.push_back("Gamma3"); 
  nodename.push_back("S035R1");
  nodetitle[inode]="G(mu- nubar(mu) nu(tau)) / G(total)";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// mu- nubar(mu) nu(tau) :: Gamma3
  ++inode;//66
  nodegammaname.push_back("Gamma150by66"); 
  nodename.push_back("S035R14");
  nodetitle[inode]="G(h- omega nu(tau)) / G(h- h- h+ pi0 nu(tau) (ex. K0))";
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- omega nu(tau) :: Gamma151
#endif
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// K- eta nu(tau) :: Gamma128
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- pi0 omega nu(tau) :: Gamma152
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
#if defined USING_NBASE31
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_den_parm[inode].push_back(803);       node_den_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//67
  nodegammaname.push_back("Gamma149"); 
  nodename.push_back("S035R15"); 
  nodetitle[inode]="G(h- omega >=0 neutrals nu(tau)) / G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- pi0 omega nu(tau) :: Gamma152
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//68
  nodegammaname.push_back("Gamma5"); 
  nodename.push_back("S035R2");
  nodetitle[inode]="G(e- nubar(e) nu(tau)) / G(total)";
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  ++inode;//69
  nodegammaname.push_back("Gamma19"); 
  nodename.push_back("S035R20");
  nodetitle[inode]="G(h- 2pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  ++inode;//70
  nodegammaname.push_back("Gamma26");
  nodename.push_back("S035R21");
  nodetitle[inode]="G(h- 3pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE1    ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  ++inode;//71
  nodegammaname.push_back("Gamma150");
  nodename.push_back("S035R23");
  nodetitle[inode]="G(h- omega nu(tau)) / G(total)";
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- omega nu(tau) :: Gamma151
#endif
  ++inode;//72
  nodegammaname.push_back("Gamma2");
  nodename.push_back("S035R24");
  nodetitle[inode]="G((particles)- >=0 neutrals >=0 K(L)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// mu- nubar(mu) nu(tau) :: Gamma3
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES1 ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PIZSQUAREPlus1]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_KS2PIZ           ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES2 ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// K- K0 nu(tau) :: Gamma37
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES2 ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 pi- eta nu(tau) :: Gamma132 
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KSKL_KS2PIZ  ]);// K- phi nu(tau) (phi->KS(->pi0pi0)KL) :: Gamma801
#endif
#endif
  ++inode;//73
  nodegammaname.push_back("Gamma31");
  nodename.push_back("S035R26");
  nodetitle[inode]="G(K- >=0 pi0 >=0 K0 >=0 gamma nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES1 ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KSKL_KS2PIZ  ]);// K- phi nu(tau) (phi->KS(->pi0pi0)KL) :: Gamma801
#endif
  ++inode;//74
  nodegammaname.push_back("Gamma32");
  nodename.push_back("S035R27");
  nodetitle[inode]="G(K- >=1 (pi0 or K0 or gamma) nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES1 ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KSKL_KS2PIZ  ]);// K- phi nu(tau) (phi->KS(->pi0pi0)KL) :: Gamma801
#endif
  ++inode;//75
  nodegammaname.push_back("Gamma56");
  nodename.push_back("S035R28");
  nodetitle[inode]="G(h- h- h+ nu(tau)) / G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIP         ]);// K- omega nu(tau) :: Gamma151
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
#endif
#endif
  ++inode;//76
  nodegammaname.push_back("Gamma63");
  nodename.push_back("S035R30");
  nodetitle[inode]="G(h- h- h+ >=1 neutrals nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PI_X_2PIZ_X_2 ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_KS2PI            ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE31
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OMPIMPIPPIZ      ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_Gam132_3PI3PIZ_KL]);// Kbar0 pi- eta nu(tau) -> 3pi- 3pi0, 3pi- pi0 KL0 :: Gamma132 
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KSKL_KS2PI   ]);// K- phi nu(tau) (phi->KS(->pi-pi+)KL) :: Gamma801
#endif
#endif
  ++inode;//77
  nodegammaname.push_back("Gamma54");
  nodename.push_back("S035R31");
  nodetitle[inode]="G(h- h- h+ >=0 neutrals >=0 K(L)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PI_X_2PIZ_X_2 ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_KS2PI            ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_Gam132_3PI3PIZ_KL]);// Kbar0 pi- eta nu(tau) -> 3pi- 3pi0, 3pi- pi0 KL0 :: Gamma132 
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(coef[C_KS2PIZby2PlusHalf]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KK_KS2PI     ]);// K- phi nu(tau) (phi->K-K+, KS(->pi-pi+)KL) :: Gamma801
#endif
#endif
  ++inode;//78
  nodegammaname.push_back("Gamma126");
  nodename.push_back("S035R32");
  nodetitle[inode]="G(pi- pi0 eta nu(tau)) / G(total)";
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 eta nu(tau) :: Gamma126
  ++inode;//79
  nodegammaname.push_back("Gamma102");
  nodename.push_back("S035R33");
  nodetitle[inode]="G(3h- 2h+ >=0 neutrals nu(tau) (ex. K(S)0 --> pi- pi+) (``5-prong'')) / G(total)";
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ nu(tau) (ex. K0) :: Gamma103
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ pi0 nu(tau) (ex. K0) :: Gamma104
  ++inode;//80
  nodegammaname.push_back("Gamma79");
  nodename.push_back("S035R34");
  nodetitle[inode]="G(K- h- h+ >=0 neutrals nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
#if defined USING_NBASE31 
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KK_KS2PI     ]);// K- phi nu(tau) (phi->K-K+, KS(->pi-pi+)KL) :: Gamma801
#endif
#endif
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//81
  nodegammaname.push_back("Gamma103");
  nodename.push_back("S035R38");
  nodetitle[inode]="G(3h- 2h+ nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ nu(tau) (ex. K0) :: Gamma103
  ++inode;//82
  nodegammaname.push_back("Gamma104");
  nodename.push_back("S035R39");
  nodetitle[inode]="G(3h- 2h+ pi0 nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ pi0 nu(tau) (ex. K0) :: Gamma104
  ++inode;//83
  nodegammaname.push_back("Gamma93");
  nodename.push_back("S035R40");
  nodetitle[inode]="G(pi- K- K+ nu(tau)) / G(total)";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  ++inode;//84
  nodegammaname.push_back("Gamma82");
  nodename.push_back("S035R41");
  nodetitle[inode]="G(K- pi- pi+ >=0 neutrals nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETACHARGEDMODES  ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- pi0 K0 nu(tau) :: Gamma42
#if defined USING_NBASE31 
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM3PIPlus2PI     ]);// K- omega nu(tau) :: Gamma151
#endif
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIby2         ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//85
  nodegammaname.push_back("Gamma11");
  nodename.push_back("S035R42");
  nodetitle[inode]="G(h- >=1 neutrals nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES1 ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PIZSQUARED    ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES2 ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE31
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_OM2PIZGAMMA      ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETANEUTRALMODES2 ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//86
  nodegammaname.push_back("Gamma7");
  nodename.push_back("S035R43");
  nodetitle[inode]="G(h- >=0 K(L)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_HALF             ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_HALF             ]);// K- K0 nu(tau) :: Gamma37
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  ++inode;//87
  nodegammaname.push_back("Gamma17");
  nodename.push_back("S035R44");
  nodetitle[inode]="G(h- >=2 pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE1    ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PIZSQUARED    ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE2    ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- K0 nu(tau) :: Gamma37
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE2    ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//88
  nodegammaname.push_back("Gamma37");
  nodename.push_back("S035R46");
  nodetitle[inode]="G(K- K0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
  ++inode;//89
  nodegammaname.push_back("Gamma3by5");
  nodename.push_back("S035R5");
  nodetitle[inode]="G(mu- nubar(mu) nu(tau)) / G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// mu- nubar(mu) nu(tau) :: Gamma3
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  ++inode;//90
  nodegammaname.push_back("Gamma9");
  nodename.push_back("S035R6");
  nodetitle[inode]="G(pi- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  ++inode;//91
  nodegammaname.push_back("Gamma10");
  nodename.push_back("S035R7");
  nodetitle[inode]="G(K- nu(tau)) / G(total)";
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  ++inode;//92
  nodegammaname.push_back("Gamma14");
  nodename.push_back("S035R8");
  nodetitle[inode]="G(pi- pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  ++inode;//93
  nodegammaname.push_back("Gamma13"); 
  nodename.push_back("S035R84"); 
  nodetitle[inode]="G(h- pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  ++inode;//94
  nodegammaname.push_back("Gamma24");
  nodename.push_back("S035R97");
  nodetitle[inode]="G(h- >= 3pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE1    ]);// K- eta nu(tau) :: Gamma128
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_KS2PIZby2        ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_KS2PIZSQUARED    ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE2    ]);// pi- pi0 eta nu(tau) :: Gamma126
#if defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE2    ]);// K- pi0 eta nu(tau) :: Gamma130
#endif
  ++inode;//95
  nodegammaname.push_back("Gamma9by5");
  nodename.push_back("S035Y01");
  nodetitle[inode]="G(pi- nu(tau)) / G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  ++inode;//96
  nodegammaname.push_back("Gamma10by5");
  nodename.push_back("S035Y02");
  nodetitle[inode]="G(K- nu(tau)) / G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  ++inode;//97
  nodegammaname.push_back("Gamma10by9");
  nodename.push_back("S035Y03");
  nodetitle[inode]="G(K- nu(tau)) / G(pi- nu(tau))";
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  node_den_parm[inode].push_back(12 );       node_den_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
#if defined USING_NBASE39 || defined USING_NBASE40
  ++inode;//98
  nodegammaname.push_back("Gamma130");
  nodename.push_back("S035C27");
  nodetitle[inode]="G(K- pi0 eta nu(tau)) / G(total)";
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 eta nu(tau) :: Gamma130
  ++inode;//99
  nodegammaname.push_back("Gamma132");
  nodename.push_back("S035C28");
  nodetitle[inode]="G(Kbar0 pi- eta nu(tau)) / G(total)";
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- eta nu(tau) :: Gamma132
  ++inode;//100
  nodegammaname.push_back("Gamma43");
  nodename.push_back("S035C37");
  nodetitle[inode]="G(Kbar0 pi- >=1 pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
  ++inode;//101
  nodegammaname.push_back("Gamma44");
  nodename.push_back("S035B98");
  nodetitle[inode]="G(Kbar0 pi- 2pi0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
  ++inode;//102
  nodegammaname.push_back("Gamma53");
  nodename.push_back("S035C5");
  nodetitle[inode]="G(Kbar0 h- h- h+ nu(tau)) / G(total)";
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
  ++inode;//103
  nodegammaname.push_back("Gamma800");
  nodename.push_back("S035Z01");
  nodetitle[inode]="G(pi- omega nu(tau)) / G(total)";
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- omega nu(tau) :: Gamma800
  ++inode;//104
  nodegammaname.push_back("Gamma151");
  nodename.push_back("S035C61");
  nodetitle[inode]="G(K- omega nu(tau)) / G(total)";
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- omega nu(tau) :: Gamma151
  ++inode;//105
  nodegammaname.push_back("Gamma802");
  nodename.push_back("S035Z03");
  nodetitle[inode]="G(K- pi- pi+ nu(tau) (ex. K0, omega)) / G(total)";
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  ++inode;//106
  nodegammaname.push_back("Gamma803");
  nodename.push_back("S035Z04");
  nodetitle[inode]="G(K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta)) / G(total)";
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  ++inode;//107
  nodegammaname.push_back("Gamma136");
  nodename.push_back("S035B89");
  nodetitle[inode]="G(pi- pi- pi+ eta nu(tau) (ex. K0)) / G(total)";
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(coef[C_ETAPIMPIPPIZ     ]);// 3h- 2h+ pi0 nu(tau) (ex. K0) :: Gamma104
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ETA3PIZUSAGE1    ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  ++inode;//108
  nodegammaname.push_back("Gamma115");
  nodename.push_back("S035B60");
  nodetitle[inode]="G(K- Kstar(892)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_TWOTHIRD         ]);// pi- K- K+ nu(tau) :: Gamma93
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONETHIRD         ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  ++inode;//109
  nodegammaname.push_back("Gamma49"); 
  nodename.push_back("S035C44"); 
  nodetitle[inode]="G(pi- pi0 K0 Kbar0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(278);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 K0 Kbar0 nu(tau) :: Gamma49  
  ++inode;//110
  nodegammaname.push_back("Gamma804"); 
  nodename.push_back("S035Z05"); 
  nodetitle[inode]="G(pi- K(L)0 K(L)0 nu(tau)) / G(total)";
  node_num_parm[inode].push_back(804);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(L)0 K(L)0 nu(tau) :: Gamma804  
  ++inode;//111
  nodegammaname.push_back("Gamma805"); 
  nodename.push_back("S035Z06"); 
  nodetitle[inode]="G(a1- (-> pi- gamma) nu(tau)) / G(total)";
  node_num_parm[inode].push_back(805);       node_num_coef[inode].push_back(coef[C_ONE              ]);// a1- (-> pi- gamma) nu(tau) :: Gamma805  
#if defined USING_NBASE40
  ++inode;//112
  nodegammaname.push_back("Gamma96");
  nodename.push_back("S035C9");
  nodetitle[inode]="G(K- K- K+ nu(tau)) / G(total)";
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_PHI2KMKP         ]);// K- K- K+ nu(tau) :: Gamma96
  ++inode;//113
  nodegammaname.push_back("Gamma801");
  nodename.push_back("S035Z02");
  nodetitle[inode]="G(K- phi nu(tau) (phi->KK)) / G(total)";
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- phi nu(tau) (phi->KK) :: Gamma801
#endif
#endif
  ++inode;//nnode-1
  nodegammaname.push_back("GammaAll"); // sum of all base nodes  
  nodename.push_back("S035S01");
  nodetitle[inode]="G(total)";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// mu- nubar(mu) nu(tau) :: Gamma3
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// e- nubar(e) nu(tau) :: Gamma5
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- nu(tau) :: Gamma9
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- nu(tau) :: Gamma10
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 nu(tau) :: Gamma14
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 nu(tau) :: Gamma16
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 2pi0 nu(tau) (ex. K0) :: Gamma20
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 2pi0 nu(tau) (ex. K0) :: Gamma23
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- 3pi0 nu(tau) (ex. K0) :: Gamma27
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- 3pi0 nu(tau) (ex. K0, eta) :: Gamma28
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- 4pi0 nu(tau) (ex. K0, eta) :: Gamma30
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- nu(tau) :: Gamma35
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- K0 nu(tau) :: Gamma37
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- pi0 nu(tau) :: Gamma40
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 K0 nu(tau) :: Gamma42
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(S)0 nu(tau) :: Gamma47
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(S)0 K(L)0 nu(tau) :: Gamma48
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma62
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi- pi+ pi0 nu(tau) (ex. K0, omega) :: Gamma70
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta) :: Gamma77
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- h- h+ 3pi0 nu(tau) :: Gamma78
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ nu(tau) :: Gamma93
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K- K+ pi0 nu(tau) :: Gamma94
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ pi0 nu(tau) (ex. K0) :: Gamma104
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 eta nu(tau) :: Gamma126
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- eta nu(tau) :: Gamma128
#if defined USING_NBASE31
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0) :: Gamma85
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, eta) :: Gamma89
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- omega nu(tau) :: Gamma150
#elif defined USING_NBASE39 || defined USING_NBASE40
  node_num_parm[inode].push_back(802);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ nu(tau) (ex. K0, omega) :: Gamma802
  node_num_parm[inode].push_back(803);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta) :: Gamma803
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- omega nu(tau) :: Gamma800
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- omega nu(tau) :: Gamma151
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- pi0 eta nu(tau) :: Gamma130
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- eta nu(tau) :: Gamma132
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 pi- 2pi0 nu(tau) :: Gamma44
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(coef[C_ONE              ]);// Kbar0 h- h- h+ nu(tau) :: Gamma53
  node_num_parm[inode].push_back(278);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- pi0 K0 Kbar0 nu(tau) :: Gamma49  
  node_num_parm[inode].push_back(804);       node_num_coef[inode].push_back(coef[C_ONE              ]);// pi- K(L)0 K(L)0 nu(tau) :: Gamma804  
  node_num_parm[inode].push_back(805);       node_num_coef[inode].push_back(coef[C_ONE              ]);// a1- (-> pi- gamma) nu(tau) :: Gamma805  
#if defined USING_NBASE40
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(coef[C_ONE              ]);// K- phi nu(tau) (phi->KK) :: Gamma801
#endif
#endif
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(coef[C_ONE              ]);// h- pi0 omega nu(tau) :: Gamma152
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(coef[C_ONE              ]);// 3h- 2h+ nu(tau) (ex. K0) :: Gamma103
  //
  cout << "Read definitions of "<< inode+1 << " nodes out of nnode = " << nnode << endl;
  //
  // Count number of parameters in numerator and denominator for each node
  //
  int            node_num_npar[nnode];
  int            node_den_npar[nnode];
  for (inode=0;inode<nnode;++inode){
    node_num_npar[inode]=node_num_parm[inode].size();
    node_den_npar[inode]=node_den_parm[inode].size();
  }
  //
  // Define array of a_nodename to simplify passing the strings to a function
  // 
  char** a_nodename = new char*[nnode];  
  for (inode=0;inode<nnode;++inode) {
    a_nodename[inode] = new char[10];
    a_nodename[inode] = (char*)(nodename[inode].data()); 
    //    cout << "inode = " << inode << " nodename = " << nodename[inode] << " a_nodename = " << a_nodename[inode] << endl;
  }
  // DEFINE node_is_base
  bool node_is_base[nnode];
  for (inode=0;inode<nnode;++inode) {
    node_is_base[inode] = (node_num_npar[inode]+node_den_npar[inode])==1;
#if defined USING_NBASE40
    if (inode==N_GAMMA96) node_is_base[inode] = false; // derived node [SPECIAL CASE : Gamma96 = (1/1.699387) * Gamma801]
#endif
  }
  //
  // PRINT NODE DEFINITION 
  //
  print_node_def(nnode, a_nodename, nodetitle, nodegammaname, node_is_base, 
		 ncoef, coef, coefnames,
		 node_num_npar, node_num_parm, node_num_coef, 
		 node_den_npar, node_den_parm, node_den_coef, 
		 baseparm, basegamma, basetitle);
  //
  // DEFINE NODE DEPENDENCE ON BASE
  //
  vector<int> node_parm[nnode]; // vector of parameters for each node
  vector<int> node_quan[nnode]; // vector of quantities for each node
  int first_quan[nnode];  // first quantity measured for each node
  define_nodes(nnode, node_num_parm, node_den_parm, baseparm,
	       node_parm, // output
	       node_quan, // output
	       first_quan // output
	       );
  //
  // READ MEASUREMENT FILE
  //
  int nmeas=0;
  int nmeas_noweak;
  int measnum[200], measnode[200], imeas1, imeas2;
  string measgammaname[200], measnodename[200], expname[200], author[200], year[200], meastitle[200];
  double measvalue[200], measerror[200], corrtemp;
  double** corrmat = new double*[200]; for (i=0;i<200;++i) corrmat[i] = new double[200];
  for (imeas1=0;imeas1<200;imeas1++) for (imeas2=0;imeas2<200;imeas2++) corrmat[imeas1][imeas2] = 0;
  //
#if defined USING_NBASE31 || defined USING_NBASE39
  string measfilename = Form("s035-fit-no-babar-belle%s.data",salephhcorr.data());
#elif defined USING_NBASE40
  string measfilename = Form("s035-fit-with-babar-belle%s.data",salephhcorr.data());
#endif
  ifstream ifs(measfilename.data());
  if (!ifs.good()) {
    cout << Form("Cannot open input file : %s\n", measfilename.data());
    exit(1);
  } else {
    cout << Form("Read input measurements from : %s\n", measfilename.data());
  }
  //
  while(ifs.good()) {
    if (ifs.eof()) break;
    char firstch(' ') ; ifs.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip such lines
      ifs.ignore(256,'\n') ;
    } else if (firstch=='*') {  // measurement line
      ifs >> measnum[nmeas] >> measgammaname[nmeas] >> measnodename[nmeas]
	  >> measvalue[nmeas] >> measerror[nmeas] 
	  >> expname[nmeas] >> author[nmeas] >> year[nmeas];
      //
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[nmeas]);      
      inode=it-nodename.begin();
      measnode[nmeas]=inode;
      //
      ifs.getline(buffer,256,'\n');
      string stringbuffer=string(buffer);
      meastitle[nmeas]="";
      bool first=false;
      for (string::iterator si=stringbuffer.begin();si!=stringbuffer.end();++si){
	if (13 == (int)*si) break; // ^M is 13
	if (' ' != *si) first=true;
	if (first) meastitle[nmeas]+=*si;
      }
      //      cout << measnum[nmeas] << " " << measgammaname[nmeas] << " " << measnodename[nmeas] << " " 
      //	   << measvalue[nmeas] << " " << measerror[nmeas] << " " 
      //	   << expname[nmeas] << " " << author[nmeas] << " " << year[nmeas] << " " << meastitle[nmeas] << endl;
      ++nmeas;
    } else if (firstch=='%') {  // correlation line
      ifs >> imeas1 >> imeas2 >> corrtemp;       ifs.ignore(256,'\n') ;
      if (corrmat[imeas1-1][imeas2-1] != 0) {cout << "WATCH OUT" << endl; exit(1);}
      corrmat[imeas1-1][imeas2-1] = corrmat[imeas2-1][imeas1-1] = corrtemp*1e-2;
      //
      cout << "Correlation between imeas1 = " << Form("%4d",measnum[imeas1-1])  
	   << " [Gamma=" << Form("%8s",measgammaname[imeas1-1].data()) << ", Node=" << Form("%4d",measnode[imeas1-1]) << "] "
	   << " and imeas2 = " << Form("%4d",measnum[imeas2-1]) 
	   << " [Gamma=" << Form("%8s",measgammaname[imeas2-1].data()) << ", Node=" << Form("%4d",measnode[imeas2-1]) << "] " 
	   << " is = " <<  Form("% 10.6f",corrmat[imeas1-1][imeas2-1]*100.) << " % " << endl;
      //
    }
  }
  //
  // Process Measurements
  //
  TMatrixD MeasErrorMatrix(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix(nmeas,nmeas);
  vector<int> icorrj[200]; // vector of measurements correlated with this measurement
  int ncorrij; // number of cycles of correlated measurements
  int ifirstj[200]; // first measurement in this cycle
  vector<int> veccorrij[200]; // vector of correlated measurements per cycle
  int newnode_all;
  int nodegroup_all[200]; // mapping of each measurement to groups of nodes
  int ncycle[200]; // mapping of each measurement to the cycle number
  int n_nodegroup[200]; // of newnode_all
  int ngroup[200]; // of measurements
  //
  process_measurements(nmeas, measnode, measerror, corrmat,
		       MeasErrorMatrix,  // output
		       InvMeasErrorMatrix, //output
		       icorrj, //output
		       ncorrij, //output
		       ifirstj, //output
		       veccorrij, //output
		       newnode_all, //output
		       nodegroup_all, //output
		       ncycle, // output
		       n_nodegroup, // output
		       ngroup //output
		       );
  //
  // PRINT INFORMATION ABOUT INPUT MEASUREMENTS
  //
#if defined USING_NBASE31 || defined USING_NBASE39
  FILE *measinfofile=fopen(Form("s035-fit-no-babar-belle%s.info",salephhcorr.data()),"w");  
#elif defined USING_NBASE40
  FILE *measinfofile=fopen(Form("s035-fit-with-babar-belle%s.info",salephhcorr.data()),"w");  
#endif
  print_measinfo(measinfofile,
		 nmeas, measnode, measvalue, measerror, corrmat,
		 expname, meastitle, measgammaname,
		 a_nodename, node_num_npar, node_den_npar, node_quan, node_is_base,
		 nbase, baseparm, basegamma, basenodename, basetitle,
		 newnode_all, nodegroup_all, ncycle, ncorrij, ifirstj, icorrj, veccorrij, ngroup);
  fclose(measinfofile);
  //
  // Derivative of R [node] as a linearized sum of variables P_i [base quantities]
  //
  double node_num[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double node_den[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double node_val[nnode]; // value of each node
  vector<double> node_part[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- may change
		   node_num, // output
		   node_den, // output
		   node_val, // output
		   node_part // output
		   );
  //
  //  Prepare to obtain Fit Values
  //
  double basevalue_fit[nbase];
  double baseerror_fit[nbase];
  double** basecov_fit = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit[i] = new double[nbase];
  double** basecorr_fit = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit[i] = new double[nbase];
  double chisquared;
  int weak_none[200];  for (i=0;i<nmeas;++i) weak_none[i]=0; // all measurements are treated to be non-weak here
  //
  // COMBINE 
  //
  combine(uconstrain,
	  nmeas, weak_none, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix,
	  node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit, // output
	  baseerror_fit, // output
	  basecov_fit,   // output
	  basecorr_fit,  // output
	  chisquared     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_temp = chisquared;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    // COMBINE 
    combine(uconstrain,
	    nmeas, weak_none, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix,
	    node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit, // output
	    baseerror_fit, // output
	    basecov_fit,   // output
	    basecorr_fit,  // output
	    chisquared     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared - chisquared_temp) < 1e-3) break;
    chisquared_temp = chisquared;
  }
  //
  cout << endl << "Results from original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE             SEED      FITVAL         FITERR    TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),baseseed[ibase],
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue[nnode];     // value of each node
  vector<double> NodeValue_part[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num, // output
		   NodeValue_den, // output
		   NodeValue,     // output
		   NodeValue_part // output
		   );
  //
  double NodeError[nnode];
  double** NodeErrorMatrix = new double*[nnode]; for (i=0;i<nnode;++i) NodeErrorMatrix[i]= new double[nnode];
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part[jnode].at(jpar);
	  NodeErrorMatrix[inode][jnode] += (iquan>nbase || jquan>nbase) ? 0 :
	    ipartial*jpartial*baseerror_fit[iquan-1]*baseerror_fit[jquan-1]*basecorr_fit[iquan-1][jquan-1];
	}
      }
    }
    if (NodeErrorMatrix[inode][inode]<0) NodeErrorMatrix[inode][inode] = 0;
    //
    NodeError[inode] = TMath::Sqrt(TMath::Max(0.,NodeErrorMatrix[inode][inode]));
    //    cout << "inode = " << inode << " NodeError = " << NodeError[inode] << " NodeErrorMatrix = " << NodeErrorMatrix[inode][inode] << endl;
  }
  //
  int   nchisquare_tot=0; 
  double chisquare_tot=0; 
  double chisquare_meas[200]; // chi2 per mesurement
  int   nchisquare_meas[200];  for (im=0;im<nmeas;++im) nchisquare_meas[im]=0;// nodes per mesurement     
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_none[jm]==weakcompare) continue;
	  chisquare_meas[im]+= 
	    (measvalue[im]-NodeValue[measnode[im]])*InvMeasErrorMatrix[im][jm]*(measvalue[jm]-NodeValue[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot+=1;
	chisquare_tot+=chisquare_node;
	nchisquare_meas[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror[im], 
		  NodeValue[measnode[im]], NodeError[measnode[im]],
		  chisquare_meas[im], chisquare_node, chisquare_tot,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g \n\n", 
	       chisquare_tot, nchisquare_tot, nbase_u, nchisquare_tot - nbase_u, TMath::Prob(chisquare_tot, nchisquare_tot - nbase_u));
  //
  cout << Form("%s = (%10.6f +- %10.6f)%% (%10.6f sigma); %s = (%10.6f +- %10.6f)%% (%10.6f sigma)\n\n",
	       nodetitle[N_GAMMAALL].data(),
	       NodeValue[N_GAMMAALL]*100.,
	       uconstrain?0:NodeError[N_GAMMAALL]*100.,
	       uconstrain?0:(1-NodeValue[N_GAMMAALL])/NodeError[N_GAMMAALL],
	       nodetitle[N_GAMMA110].data(),
	       NodeValue[N_GAMMA110]*100.,
	       NodeError[N_GAMMA110]*100.,
	       NodeValue[N_GAMMA110]/NodeError[N_GAMMA110]);
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile[2];
  for (int p=0;p<2;++p){
    if (p==0) measfile[0]=fopen(Form("combos_measurements_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile[1]=fopen(Form("alucomb_measurements_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_none, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror, corrmat,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile[2];
  for (int p=0;p<2;++p){
    if (p==0) avefile[0]=fopen(Form("combos_average_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile[1]=fopen(Form("alucomb_average_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile[p], "INCLUDE combos_measurements_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile[p], "INCLUDE alucomb_measurements_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_none, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, node_is_base,
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  //  Prepare to obtain Fit Values for non-weak measurements
  //
  double basevalue_fit_noweak[nbase];
  double baseerror_fit_noweak[nbase];
  double** basecov_fit_noweak = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_noweak[i] = new double[nbase];
  double** basecorr_fit_noweak = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_noweak[i] = new double[nbase];
  double chisquared_noweak;
  //
  // Define weak nodes
  //
  double nsig[200];
  int weak[200];
  for (i=0;i<nmeas;++i){
    if (nchisquare_meas[i] > 0 ) {
      nsig[i] = measerror[i]/((sqrt(nchisquare_meas[i]))*NodeError[measnode[i]]);
    } else {
      nsig[i] = 999999;
    }
    weak[i] = (icorrj[i].size()==1) && (nsig[i] > 3.); // only for uncorrelated measurements
  }
  //
  // COMBINE 
  //
  combine(uconstrain,
	  nmeas, weak, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix,
	  node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_noweak, // output
	  baseerror_fit_noweak, // output
	  basecov_fit_noweak,   // output
	  basecorr_fit_noweak,  // output
	  chisquared_noweak     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_noweak_temp = chisquared_noweak;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_noweak[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE 
    //
    combine(uconstrain,
	    nmeas, weak, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix,
	    node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_noweak, // output
	    baseerror_fit_noweak, // output
	    basecov_fit_noweak,   // output
	    basecorr_fit_noweak,  // output
	    chisquared_noweak     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared_noweak, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_noweak - chisquared_noweak_temp) < 1e-3) break;
    chisquared_noweak_temp = chisquared_noweak;
  }
  //
  cout << endl << "Comparison of Results from fit with [non-weak measurements only] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE      ORIG_FITVAL    ORIG_FITERR    SCAL_FITVAL    SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %14.7f %6.3f  %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_noweak[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_noweak[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_noweak[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_noweak[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_noweak[nnode];     // value of each node
  vector<double> NodeValue_part_noweak[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_noweak, // output
		   NodeValue_den_noweak, // output
		   NodeValue_noweak,     // output
		   NodeValue_part_noweak // output
		   );
  //
  double NodeError_noweak[nnode];
  double** NodeErrorMatrix_noweak = new double*[nnode]; for (i=0;i<nnode;++i) NodeErrorMatrix_noweak[i]= new double[nnode];
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_noweak[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_noweak[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_noweak[jnode].at(jpar);
	  NodeErrorMatrix_noweak[inode][jnode] += (iquan>nbase || jquan>nbase) ? 0 :
	    ipartial*jpartial*baseerror_fit_noweak[iquan-1]*baseerror_fit_noweak[jquan-1]*basecorr_fit_noweak[iquan-1][jquan-1];
	}
      }
    }
    if (NodeErrorMatrix_noweak[inode][inode]<0) NodeErrorMatrix_noweak[inode][inode] = 0;
    //
    NodeError_noweak[inode] = TMath::Sqrt(TMath::Max(0.,NodeErrorMatrix_noweak[inode][inode]));
  }
  //
  int   nchisquare_tot_noweak=0; 
  double chisquare_tot_noweak=0; 
  double chisquare_meas_noweak[200]; // chi2 per mesurement
  int   nchisquare_meas_noweak[200]; for (im=0;im<nmeas;++im) nchisquare_meas_noweak[im]=0;// nodes per mesurement
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_noweak[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak[jm]==weakcompare) continue;
	  chisquare_meas_noweak[im]+= 
	    (measvalue[im]-NodeValue_noweak[measnode[im]])*InvMeasErrorMatrix[im][jm]*(measvalue[jm]-NodeValue_noweak[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_noweak[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_noweak+=1;
	chisquare_tot_noweak+=chisquare_node;
	nchisquare_meas_noweak[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror[im], 
		  NodeValue_noweak[measnode[im]], NodeError_noweak[measnode[im]],
		  chisquare_meas_noweak[im], chisquare_node, chisquare_tot_noweak,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g \n\n", 
	       chisquare_tot_noweak, nchisquare_tot_noweak, 
	       nbase_u, nchisquare_tot_noweak - nbase_u, TMath::Prob(chisquare_tot_noweak, nchisquare_tot_noweak - nbase_u));
  //
  cout << Form("%s = (%10.6f +- %10.6f)%% (%10.6f sigma); %s = (%10.6f +- %10.6f)%% (%10.6f sigma)\n\n",
	       nodetitle[N_GAMMAALL].data(),
	       NodeValue_noweak[N_GAMMAALL]*100.,
	       uconstrain?0:NodeError_noweak[N_GAMMAALL]*100.,
	       uconstrain?0:(1-NodeValue_noweak[N_GAMMAALL])/NodeError_noweak[N_GAMMAALL],
	       nodetitle[N_GAMMA110].data(),
	       NodeValue_noweak[N_GAMMA110]*100.,
	       NodeError_noweak[N_GAMMA110]*100.,
	       NodeValue_noweak[N_GAMMA110]/NodeError_noweak[N_GAMMA110]);
  //
  cout << "Summary of fit with [non-weak measurements only]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE             SEED      FITVAL         FITERR     SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %6.3f %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),baseseed[ibase],
		 basevalue_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // PDG-style Scale Factor Calculation
  //
  printf ("ncorrij = %d \n",ncorrij);
  for (i=0;i<ncorrij;++i){
    printf ("i = %d ifirstj[i]+1 = %3d ",i,ifirstj[i]+1);
    printf ("veccorrij[i].size() = %2d veccorrij[i][j]+1 = ",veccorrij[i].size());
    for (j=0;j<veccorrij[i].size();++j) printf ("%d ",veccorrij[i][j]+1);
    printf ("GammaName = "); 
    for (j=0;j<veccorrij[i].size();++j) printf ("%s ",measgammaname[veccorrij[i][j]].data());
    printf ("\n");
  }
  //
  vector<bool> VectorOfEigenVal_AllPositive;
  vector<double> VectorOfPullMag;
  vector<TMatrixD> VectorOfMeasScaledErrorMatrix;
  //
  vector<bool> VectorOfEigenVal_AllPositive_noweak;
  vector<double> VectorOfPullMag_noweak;
  vector<TMatrixD> VectorOfMeasScaledErrorMatrix_noweak;
  //
  for (int iweak=0;iweak<2;++iweak) {
    for (int ij=0;ij<ncorrij;++ij) {
      int nij=veccorrij[ij].size();
      TMatrixD ThisFitErrorMatrix(nij,nij);
      TMatrixD ThisMeasErrorMatrix(nij,nij);
      for (i=0;i<nij;++i){
	im=veccorrij[ij][i];
	inode=measnode[im];
	for (j=0;j<nij;++j){
	  jm=veccorrij[ij][j];
	  jnode=measnode[jm];
	  ThisMeasErrorMatrix[i][j] = (i==j) ? measerror[im]*measerror[im] : corrmat[im][jm] * measerror[im] * measerror[jm];
	  ThisFitErrorMatrix[i][j]  = (iweak==0) ? NodeErrorMatrix[inode][jnode] : NodeErrorMatrix_noweak[inode][jnode];
	  //	  if (iweak==1&&(ij==1||ij==7)) {
	  //	    cout << i << " " << j << " " << inode << " " << jnode << " " << nodegammaname[inode] << " " << nodegammaname[jnode] << " "
	  //		 << "NodeErrorMatrix: " << NodeErrorMatrix[inode][jnode] << " " << NodeErrorMatrix_noweak[inode][jnode] << " "
	  //		 << NodeErrorMatrix[inode][jnode]/NodeErrorMatrix_noweak[inode][jnode] << " "
	  //		 << endl;
	}
      }
      //
      // CMatrix = (MeasErrorMatrix - FitErrorMatrix)^{-1}
      //
      TMatrixD InvCMatrix(nij,nij);
      for (i=0;i<nij;++i) {
	for (j=0;j<nij;++j) {
	  InvCMatrix[i][j] = ThisMeasErrorMatrix[i][j] - ThisFitErrorMatrix[i][j];
	}
      }
      TMatrixD CMatrix = InvCMatrix;
      Double_t c_det;
      CMatrix.Invert(&c_det);
      //
      TVectorD EigenVal;
      TMatrixD EigenVecT = CMatrix.EigenVectors(EigenVal);
      TMatrixD EigenVec(nij,nij); EigenVec.Transpose(EigenVecT);
      cout << "EigenVal : "; for (i=0;i<nij;++i) cout << EigenVal[i] << " ";  cout << endl;
      bool EigenVal_AllPositive=true; for (i=0;i<nij;++i) if (EigenVal[i]<0) EigenVal_AllPositive=false;
      //
      double PullMag=0;
      double PullMag2=0;
      double PullVector[nij];
      for (i=0;i<nij;++i) {
	PullVector[i] = 0;
      }
      TMatrixD ThisMeasScaledErrorMatrix(nij,nij);
      //
      if (EigenVal_AllPositive) {
	TMatrixD DiagonalMatrix(nij,nij);
	TMatrixD DiagonalMatrix_sqrt(nij,nij);
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    if (i==j) {
	      DiagonalMatrix(i,j) = EigenVal[i];
	      DiagonalMatrix_sqrt(i,j) = TMath::Sqrt(EigenVal[i]);
	    } else {
	      DiagonalMatrix(i,j) = 0;
	      DiagonalMatrix_sqrt(i,j) = 0;
	    }
	  }
	}
	//
	TMatrixD QMatrix = EigenVecT * DiagonalMatrix_sqrt * EigenVec; // Q = sqrt(C)
	//
	double FMinusM[nij];
	for (i=0;i<nij;++i) {
	  im=veccorrij[ij][i];
	  FMinusM[i] = (iweak==0) ? NodeValue[measnode[im]]-measvalue[im] : NodeValue_noweak[measnode[im]]-measvalue[im]; 
	}
	//
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    PullVector[i] += QMatrix[i][j] * FMinusM[j];
	  }
	  PullMag2 += PullVector[i]*PullVector[i];
	}
	PullMag = TMath::Sqrt(PullMag2);
	//
	double UnitPullVector[nij];
	for (i=0;i<nij;++i) {
	  UnitPullVector[i] = (PullMag>0)? PullVector[i]/PullMag : 0;
	}
	TMatrixD ThisMeasCorrMatrix(nij,nij);
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    if (i==j) {
	      ThisMeasCorrMatrix[i][j] = 1; // by construction diagonal elements of corrmat has not been filled
	    } else {
	      im=veccorrij[ij][i];
	      jm=veccorrij[ij][j];
	      ThisMeasCorrMatrix[i][j]=corrmat[im][jm];
	    }
	  }
	}
	//
	double CM=0;
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    CM += UnitPullVector[i]*ThisMeasCorrMatrix[i][j]*UnitPullVector[j];
	  }
	}
	double SFactor = (PullMag>1) ? (PullMag2 - 1) * CM : 0; 
	//
	for (i=0;i<nij;++i) {
	  im=veccorrij[ij][i];
	  for (j=0;j<nij;++j) {
	    jm=veccorrij[ij][j];
	    ThisMeasScaledErrorMatrix[i][j] = ThisMeasErrorMatrix[i][j] + UnitPullVector[i] * measerror[im] * UnitPullVector[j] * measerror[jm] * SFactor;
	  }
	}
      }
      //
      cout << "iweak = " << iweak << " ij = " << ij << " PullMag = " << PullMag << " PullVector : "; for (i=0;i<nij;++i) cout << PullVector[i] << " "; cout << endl;
      //
      if (iweak==0) {
	VectorOfEigenVal_AllPositive.push_back(EigenVal_AllPositive);
	if (ij==ncorrij-1) {
	  bool status = true;
	  for (i=0;i<ncorrij;++i){
	    if (!VectorOfEigenVal_AllPositive.at(i)) status = false;
	  }
	  if (status) {
	    cout << "iweak = " << iweak << " VectorOfEigenVal_AllPositive has status = TRUE " << endl << endl;
	  } else {
	    cout << "iweak = " << iweak << " VectorOfEigenVal_AllPositive has status = FALSE" << endl << endl;
	  }
	}
	VectorOfPullMag.push_back(PullMag);
	VectorOfMeasScaledErrorMatrix.push_back(ThisMeasScaledErrorMatrix);
      } else {
	VectorOfEigenVal_AllPositive_noweak.push_back(EigenVal_AllPositive);
	if (ij==ncorrij-1) {
	  bool status = true;
	  for (i=0;i<ncorrij;++i){
	    if (!VectorOfEigenVal_AllPositive_noweak.at(i)) status = false;
	  }
	  if (status) {
	    cout << "iweak = " << iweak << " VectorOfEigenVal_AllPositive has status = TRUE " << endl << endl;
	  } else {
	    cout << "iweak = " << iweak << " VectorOfEigenVal_AllPositive has status = FALSE" << endl << endl;
	  }
	}
	VectorOfPullMag_noweak.push_back(PullMag);
	VectorOfMeasScaledErrorMatrix_noweak.push_back(ThisMeasScaledErrorMatrix);
      }
      //
    } //end of loop over ij
  } // end of loop over iweak
  //
  // With all measurements
  //
  double pullav_meas[200]; for (im=0;im<nmeas;++im) pullav_meas[im]=0; // pullav per measurement
  for (inode=0;inode<=newnode_all;++inode){
    int   npullsq_node=0;
    double pullsq_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	double pullsquare=0;
	if ( icorrj[im].size() > 1 ) { // first find which cycle this correlated measurement belongs to, then just fetch the result
	  pullsquare = TMath::Power(VectorOfPullMag[ncycle[im]],2);
	} else {
	  bool singular_measurement = TMath::Abs(measvalue[im] - NodeValue[measnode[im]]) < 1e-6 && TMath::Abs(measerror[im] - NodeError[measnode[im]]) < 1e-6;
	  pullsquare = singular_measurement ? 0 : TMath::Power((measvalue[im]-NodeValue[measnode[im]]),2)/
	    (TMath::Power(measerror[im],2) - TMath::Power(NodeError[measnode[im]],2));
	}
	pullsq_node+=pullsquare;
	npullsq_node+=1;
      }
    }
    if (npullsq_node>0) pullsq_node/=npullsq_node;
    //
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	pullav_meas[im]=TMath::Sqrt(pullsq_node);
	//cout << "inode = " << inode << " im+1 = " << im+1 << " npullsq_node = " << npullsq_node << " pullav_meas = " << pullav_meas[im] << endl;
      }
    }
  }
  //
  double measerror_scaled[200];
  for (i=0;i<nmeas;++i){
    if (icorrj[i].size() > 1 && VectorOfPullMag[ncycle[i]]>1) {
      double temperr2=measerror[i]*measerror[i]; // over-written below
      for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	if (i == veccorrij[ncycle[i]][itemp]) {
	  temperr2=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(itemp,itemp);
	  break;
	}
      }
      measerror_scaled[i] = TMath::Sqrt(temperr2);
    } else if (icorrj[i].size()==1 && pullav_meas[i]>1 ) { 
      measerror_scaled[i] = measerror[i]*pullav_meas[i];
    } else {
      measerror_scaled[i] = measerror[i];
    }
  }
  //
  double** corrmat_scaled = new double*[200]; for (i=0;i<200;++i) corrmat_scaled[i] = new double[200];
  for (int i=0;i<nmeas;++i){
    for (int j=0;j<nmeas;++j) {
      double tempcorrij=corrmat[i][j]; // over-written below
      if (icorrj[i].size() > 1 && VectorOfPullMag[ncycle[i]]>1) {
	for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	  if (i == veccorrij[ncycle[i]][itemp]) {
	    for (int jtemp=0; jtemp < veccorrij[ncycle[i]].size(); ++jtemp) {
	      if (j == veccorrij[ncycle[i]][jtemp]) {
		double temperr_i=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(itemp,itemp);
		double temperr_j=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(jtemp,jtemp);
		double temperr_ij=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(itemp,jtemp);
		tempcorrij = temperr_ij / TMath::Sqrt(temperr_i * temperr_j);
		break;
	      }
	    }
	    break;
	  }
	}
      }
      corrmat_scaled[i][j] = i==j ? 0 : tempcorrij;
      // cout << " i+1 = " << i+1 << " j+1 = " << j+1 
      //   << " icorrij = " << icorrj[i].size() << " pullav_meas = " << pullav_meas[i] << " PullMag = " << VectorOfPullMag[ncycle[i]] 
      //   << " corrmat = " << corrmat[i][j] << " corrmat_scaled = " << corrmat_scaled[i][j]
      //   << " measerror^2 = " << measerror[i] * measerror[i] 
      //   << " measerror_scaled^2 = " << measerror_scaled[i]*measerror_scaled[i]
      //   << endl;
    }
  }
  //
  TMatrixD MeasErrorMatrix_scaled(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix_scaled(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) { // by construction diagonal elements of corrmat has not been filled
	MeasErrorMatrix_scaled(i,j)=measerror_scaled[i]*measerror_scaled[i];
      } else {
	MeasErrorMatrix_scaled(i,j)=corrmat_scaled[i][j]*measerror_scaled[i]*measerror_scaled[j]; 
      }
    }
  }
  InvMeasErrorMatrix_scaled = MeasErrorMatrix_scaled;
  Double_t determinant_scaled;
  InvMeasErrorMatrix_scaled.Invert(&determinant_scaled);  
  //  cout << Form("Determinant of MeasErrorMatrix_scaled = %10.5g\n", determinant_scaled);
  //
  // With non-weak measurements
  //
  double pullav_meas_noweak[200]; for (im=0;im<nmeas;++im) pullav_meas_noweak[im]=0; // pullav per measurement
  for (inode=0;inode<=newnode_all;++inode){
    int   npullsq_node_noweak=0;
    double pullsq_node_noweak=0;
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	double pullsquare=0;
	if ( icorrj[im].size() > 1 ) { // first find which cycle this correlated measurement belongs to, then just fetch the result
	  pullsquare = TMath::Power(VectorOfPullMag_noweak[ncycle[im]],2);
	} else {
	  bool singular_measurement = TMath::Abs(measvalue[im] - NodeValue[measnode[im]]) < 1e-6 && TMath::Abs(measerror[im] - NodeError[measnode[im]]) < 1e-6;
	  pullsquare = singular_measurement ? 0 : TMath::Power((measvalue[im]-NodeValue_noweak[measnode[im]]),2)/
	    (TMath::Power(measerror[im],2) - TMath::Power(NodeError_noweak[measnode[im]],2));
	}
	pullsq_node_noweak+=pullsquare;
	npullsq_node_noweak+=1;
      }
    }
    if (npullsq_node_noweak>0) pullsq_node_noweak/=npullsq_node_noweak;
    //
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	pullav_meas_noweak[im]=TMath::Sqrt(pullsq_node_noweak);
	//cout << "inode = " << inode << " im+1 = " << im+1 << " npullsq_node_noweak = " << npullsq_node_noweak << " pullav_meas_noweak = " << pullav_meas_noweak[im] << endl;
      }
    }
  }
  //
  double measerror_scaled_noweak[200];
  for (i=0;i<nmeas;++i){
    if (icorrj[i].size() > 1 && VectorOfPullMag_noweak[ncycle[i]]>1) {
      double temperr2=measerror[i]*measerror[i]; // over-written below
      for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	if (i == veccorrij[ncycle[i]][itemp]) {
	  temperr2=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(itemp,itemp);
	  break;
	}
      }
      measerror_scaled_noweak[i] = TMath::Sqrt(temperr2);
    } else if (icorrj[i].size()==1 && pullav_meas_noweak[i]>1 ) { 
      measerror_scaled_noweak[i] = measerror[i]*pullav_meas_noweak[i];
    } else {
      measerror_scaled_noweak[i] = measerror[i];
    }
  }
  //
  double** corrmat_scaled_noweak = new double*[200]; for (i=0;i<200;++i) corrmat_scaled_noweak[i] = new double[200];
  for (int i=0;i<nmeas;++i){
    for (int j=0;j<nmeas;++j) {
      double tempcorrij=corrmat[i][j]; // over-written below
      if (icorrj[i].size() > 1 && VectorOfPullMag_noweak[ncycle[i]]>1) {
	for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	  if (i == veccorrij[ncycle[i]][itemp]) {
	    for (int jtemp=0; jtemp < veccorrij[ncycle[i]].size(); ++jtemp) {
	      if (j == veccorrij[ncycle[i]][jtemp]) {
		double temperr_i=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(itemp,itemp);
		double temperr_j=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(jtemp,jtemp);
		double temperr_ij=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(itemp,jtemp);
		tempcorrij = temperr_ij / TMath::Sqrt(temperr_i * temperr_j);
		break;
	      }
	    }
	    break;
	  }
	}
      }
      corrmat_scaled_noweak[i][j] = i==j ? 0 : tempcorrij;
      // cout << " i+1 = " << i+1 << " j+1 = " << j+1 
      //   << " icorrij = " << icorrj[i].size() << " pullav_meas_noweak = " << pullav_meas_noweak[i] << " PullMag = " << VectorOfPullMag[ncycle[i]] 
      //   << " corrmat = " << corrmat[i][j] << " corrmat_scaled_noweak = " << corrmat_scaled_noweak[i][j]
      //   << " measerror^2 = " << measerror[i] * measerror[i] 
      //   << " measerror_scaled_noweak^2 = " << measerror_scaled_noweak[i]*measerror_scaled_noweak[i]
      //   << endl;
    }
  }
  //
  TMatrixD MeasErrorMatrix_scaled_noweak(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix_scaled_noweak(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) { // by construction diagonal elements of corrmat has not been filled
	MeasErrorMatrix_scaled_noweak(i,j)=measerror_scaled_noweak[i]*measerror_scaled_noweak[i];
      } else {
	MeasErrorMatrix_scaled_noweak(i,j)=corrmat_scaled_noweak[i][j]*measerror_scaled_noweak[i]*measerror_scaled_noweak[j]; 
      }
    }
  }
  InvMeasErrorMatrix_scaled_noweak = MeasErrorMatrix_scaled_noweak;
  Double_t determinant_scaled_noweak;
  InvMeasErrorMatrix_scaled_noweak.Invert(&determinant_scaled_noweak);  
  //  cout << Form("Determinant of MeasErrorMatrix_scaled = %10.5g\n", determinant_scaled_noweak);
  //
  //
  //  Prepare to obtain Fit Values for non-weak measurements [scaled]
  //
  double basevalue_fit_noweak_scaled[nbase];
  double baseerror_fit_noweak_scaled[nbase];
  double** basecov_fit_noweak_scaled = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_noweak_scaled[i] = new double[nbase];
  double** basecorr_fit_noweak_scaled = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_noweak_scaled[i] = new double[nbase];
  double chisquared_noweak_scaled;
  //
  // Define weak nodes
  //
  double nsig_scaled[200];
  int weak_scaled[200];
  for (i=0;i<nmeas;++i){
    if (nchisquare_meas_noweak[i] > 0 ) {
      nsig_scaled[i] = measerror[i]/(TMath::Sqrt(nchisquare_meas_noweak[i])*NodeError_noweak[measnode[i]]);
    } else {
      nsig_scaled[i] = 999999;
    }
    weak_scaled[i] = (icorrj[i].size()==1) && (nsig_scaled[i] > 3.); // only for uncorrelated measurements
    //    weak_scaled[i] = weak[i];
    //
    //    cout << i 
    //	 << " nchisquare_meas = " << nchisquare_meas[i] << " nchisquare_meas_noweak = " << nchisquare_meas_noweak[i] 
    //	 << " nsig = " << nsig[i] << " nsig_scaled = " << nsig_scaled[i] 
    //	 << " weak = " << weak[i] << " weak_scaled = " << weak_scaled[i] 
    //	 << endl;
  }
  //
  // COMBINE
  //
  combine(uconstrain,
	  nmeas, weak_scaled, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix_scaled_noweak,
	  node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_noweak_scaled, // output
	  baseerror_fit_noweak_scaled, // output
	  basecov_fit_noweak_scaled,   // output
	  basecorr_fit_noweak_scaled,  // output
	  chisquared_noweak_scaled     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak_scaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_scaled,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_noweak_scaled_temp = chisquared_noweak_scaled;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_noweak_scaled[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE
    //
    combine(uconstrain,
	    nmeas, weak_scaled, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix_scaled_noweak,
	    node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_noweak_scaled, // output
	    baseerror_fit_noweak_scaled, // output
	    basecov_fit_noweak_scaled,   // output
	    basecorr_fit_noweak_scaled,  // output
	    chisquared_noweak_scaled     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak_scaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_scaled,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_noweak_scaled - chisquared_noweak_scaled_temp) < 1e-3) break;
    chisquared_noweak_scaled_temp = chisquared_noweak_scaled;
  }
  //
  cout << endl << "Comparison of Results from fit with [non-weak measurements only] and [errors inflated with PDG-style scale factors] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE      ORIG_FITVAL    ORIG_FITERR    SCAL_FITVAL    SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %14.7f %6.3f  %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_noweak_scaled[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_noweak_scaled[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_noweak_scaled[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_noweak_scaled[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_noweak_scaled[nnode];     // value of each node
  vector<double> NodeValue_part_noweak_scaled[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_noweak_scaled, // output
		   NodeValue_den_noweak_scaled, // output
		   NodeValue_noweak_scaled,     // output
		   NodeValue_part_noweak_scaled // output
		   );
  //
  double NodeError_noweak_scaled[nnode];
  double** NodeErrorMatrix_noweak_scaled = new double*[nnode]; for (i=0;i<nnode;++i) NodeErrorMatrix_noweak_scaled[i]= new double[nnode];
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_noweak_scaled[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_noweak_scaled[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_noweak_scaled[jnode].at(jpar);
	  NodeErrorMatrix_noweak_scaled[inode][jnode] += (iquan>nbase || jquan>nbase) ? 0 :
	    ipartial*jpartial*baseerror_fit_noweak_scaled[iquan-1]*baseerror_fit_noweak_scaled[jquan-1]*basecorr_fit_noweak_scaled[iquan-1][jquan-1];
	}
      }
    }
    if (NodeErrorMatrix_noweak_scaled[inode][inode]<0) NodeErrorMatrix_noweak_scaled[inode][inode] = 0;
    //
    NodeError_noweak_scaled[inode] = TMath::Sqrt(TMath::Max(0.,NodeErrorMatrix_noweak_scaled[inode][inode]));
  }
  //
  int   nchisquare_tot_noweak_scaled=0; 
  double chisquare_tot_noweak_scaled=0; 
  double chisquare_meas_noweak_scaled[200]; // chi2 per mesurement
  int   nchisquare_meas_noweak_scaled[200]; for (im=0;im<nmeas;++im) nchisquare_meas_noweak_scaled[im]=0;// nodes per mesurement
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_scaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_noweak_scaled[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_scaled[jm]==weakcompare) continue;
	  chisquare_meas_noweak_scaled[im]+= 
	    (measvalue[im]-NodeValue_noweak_scaled[measnode[im]])*InvMeasErrorMatrix_scaled_noweak[im][jm]*(measvalue[jm]-NodeValue_noweak_scaled[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_noweak_scaled[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_scaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_noweak_scaled+=1;
	chisquare_tot_noweak_scaled+=chisquare_node;
	nchisquare_meas_noweak_scaled[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror_scaled_noweak[im], 
		  NodeValue_noweak_scaled[measnode[im]], NodeError_noweak_scaled[measnode[im]],
		  chisquare_meas_noweak_scaled[im], chisquare_node, chisquare_tot_noweak_scaled,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g \n\n", 
	       chisquare_tot_noweak_scaled, nchisquare_tot_noweak_scaled, 
	       nbase_u, nchisquare_tot_noweak_scaled - nbase_u, TMath::Prob(chisquare_tot_noweak_scaled, nchisquare_tot_noweak_scaled - nbase_u));
  //
  cout << Form("%s = (%10.6f +- %10.6f)%% (%10.6f sigma); %s = (%10.6f +- %10.6f)%% (%10.6f sigma)\n\n",
	       nodetitle[N_GAMMAALL].data(),
	       NodeValue_noweak_scaled[N_GAMMAALL]*100.,
	       uconstrain?0:NodeError_noweak_scaled[N_GAMMAALL]*100.,
	       uconstrain?0:(1-NodeValue_noweak_scaled[N_GAMMAALL])/NodeError_noweak_scaled[N_GAMMAALL],
	       nodetitle[N_GAMMA110].data(),
	       NodeValue_noweak_scaled[N_GAMMA110]*100.,
	       NodeError_noweak_scaled[N_GAMMA110]*100.,
	       NodeValue_noweak_scaled[N_GAMMA110]/NodeError_noweak_scaled[N_GAMMA110]);
  //
  cout << "Summary of fit with [non-weak measurements only] and [errors inflated with PDG-style scale factors]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE             SEED      FITVAL         FITERR     SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %6.3f %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),baseseed[ibase],
		 basevalue_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile_scaled[2];
  for (int p=2;p<2;++p){
    if (p==0) measfile_scaled[0]=fopen(Form("combos_measurements_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile_scaled[1]=fopen(Form("alucomb_measurements_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile_scaled[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_scaled, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror_scaled_noweak, corrmat_scaled_noweak,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile_scaled[2];
  for (int p=2;p<2;++p){
    if (p==0) avefile_scaled[0]=fopen(Form("combos_average_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile_scaled[1]=fopen(Form("alucomb_average_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile_scaled[p], "INCLUDE combos_measurements_scaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile_scaled[p], "INCLUDE alucomb_measurements_scaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile_scaled[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_scaled, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, node_is_base,
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  // Applying Ad-Hoc Scale Factors
  //
  double measerror_rescaled[200];
  for (i=0;i<nmeas;++i){
    measerror_rescaled[i] = measerror[i];
#if defined USING_NBASE31 
    if (measnode[i]==N_GAMMA85) {// NAME = S035C21 GAMMA = 85 TITLE = G(K- pi+ pi- nu(tau) (ex. K0)) / G(total)
      double rescale = 1;//1.59;
      measerror_rescaled[i] = measerror[i] * rescale;
      cout << i << " " << measnode[i]  << " " << measgammaname[i] << " " << measnodename[i] << " " << expname[i] << " " << meastitle[i] << " " 
	   << NodeValue[measnode[i]] << " " << NodeError[measnode[i]] << " " << measvalue[i] << " " << measerror[i]*rescale << " inflated by Scale Factor = " << rescale << endl;
    }
#elif defined USING_NBASE40
    if (measnode[i]==N_GAMMA96) {// NAME = S035C9 GAMMA = 96 TITLE = G(K- K+ K- nu(tau)) / G(total)
      double rescale = 5.435276;	
      measerror_rescaled[i] = measerror[i] * rescale;
      cout << i << " " << measnode[i]  << " " << measgammaname[i] << " " << measnodename[i] << " " << expname[i] << " " << meastitle[i] << " " 
	   << NodeValue[measnode[i]] << " " << NodeError[measnode[i]] << " " << measvalue[i] << " " << measerror[i]*rescale << " inflated by Scale Factor = " << rescale << endl;
    }
#endif
  }
  //
  TMatrixD MeasErrorMatrix_rescaled(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix_rescaled(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) {
	MeasErrorMatrix_rescaled(i,j)=measerror_rescaled[i]*measerror_rescaled[i];
      } else {
	MeasErrorMatrix_rescaled(i,j)=corrmat[i][j]*measerror_rescaled[i]*measerror_rescaled[j]; // by construction diagonal elements of corrmat has not been filled
      }
    }
  }
  InvMeasErrorMatrix_rescaled = MeasErrorMatrix_rescaled;
  Double_t determinant_rescaled;
  InvMeasErrorMatrix_rescaled.Invert(&determinant_rescaled);  
  //  cout << Form("Determinant of MeasErrorMatrix_rescaled = %10.5g\n", determinant_rescaled);
  //
  //  Prepare to obtain Fit Values [rescaled]
  //
  double basevalue_fit_rescaled[nbase];
  double baseerror_fit_rescaled[nbase];
  double** basecov_fit_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_rescaled[i] = new double[nbase];
  double** basecorr_fit_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_rescaled[i] = new double[nbase];
  double chisquared_rescaled;
  int weak_none_rescaled[200];  for (i=0;i<nmeas;++i) weak_none_rescaled[i]=0; // all measurements are treated to be non-weak here
  //
  // Derivative of R [node] as a linearized sum of variables P_i [base quantities]
  //
  for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=baseseed_orig[ibase]; // <-- original baseseed re-used
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- may change
		   node_num, // output
		   node_den, // output
		   node_val, // output
		   node_part // output
		   );
  //
  // COMBINE
  //
  combine(uconstrain,
	  nmeas, weak_none_rescaled, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix_rescaled,
	  node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_rescaled, // output
	  baseerror_fit_rescaled, // output
	  basecov_fit_rescaled,   // output
	  basecorr_fit_rescaled,  // output
	  chisquared_rescaled     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_rescaled,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  //
  double chisquared_rescaled_temp = chisquared_rescaled;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_rescaled[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE
    //
    combine(uconstrain,
	    nmeas, weak_none_rescaled, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix_rescaled,
	    node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_rescaled, // output
	    baseerror_fit_rescaled, // output
	    basecov_fit_rescaled,   // output
	    basecorr_fit_rescaled,  // output
	    chisquared_rescaled     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_rescaled,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_rescaled - chisquared_rescaled_temp) < 1e-3) break;
    chisquared_rescaled_temp = chisquared_rescaled;
  }
  //
  cout << endl << "Comparison of Results from fit with [errors rescaled in a Ad-Hoc style for kkk only] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE      ORIG_FITVAL    ORIG_FITERR    SCAL_FITVAL    SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %14.7f %6.3f  %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_rescaled[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_rescaled[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_rescaled[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_rescaled[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_rescaled[nnode];     // value of each node
  vector<double> NodeValue_part_rescaled[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_rescaled, // output
		   NodeValue_den_rescaled, // output
		   NodeValue_rescaled,     // output
		   NodeValue_part_rescaled // output
		   );
  //
  double NodeError_rescaled[nnode];
  double** NodeErrorMatrix_rescaled = new double*[nnode]; for (i=0;i<nnode;++i) NodeErrorMatrix_rescaled[i]= new double[nnode];
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_rescaled[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_rescaled[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_rescaled[jnode].at(jpar);
	  NodeErrorMatrix_rescaled[inode][jnode] += (iquan>nbase || jquan>nbase) ? 0 :
	    ipartial*jpartial*baseerror_fit_rescaled[iquan-1]*baseerror_fit_rescaled[jquan-1]*basecorr_fit_rescaled[iquan-1][jquan-1];
	}
      }
    }
    if (NodeErrorMatrix_rescaled[inode][inode]<0) NodeErrorMatrix_rescaled[inode][inode] = 0;
    //
    NodeError_rescaled[inode] = TMath::Sqrt(TMath::Max(0.,NodeErrorMatrix_rescaled[inode][inode]));
  }
  //
  int   nchisquare_tot_rescaled=0; 
  double chisquare_tot_rescaled=0; 
  double chisquare_meas_rescaled[200]; // chi2 per mesurement
  int   nchisquare_meas_rescaled[200];  for (im=0;im<nmeas;++im) nchisquare_meas_rescaled[im]=0;// nodes per mesurement     
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_none_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_rescaled[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_none_rescaled[jm]==weakcompare) continue;
	  chisquare_meas_rescaled[im]+= 
	    (measvalue[im]-NodeValue_rescaled[measnode[im]])*InvMeasErrorMatrix_rescaled[im][jm]*(measvalue[jm]-NodeValue_rescaled[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_rescaled[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_none_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_rescaled+=1;
	chisquare_tot_rescaled+=chisquare_node;
	nchisquare_meas_rescaled[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror_rescaled[im], 
		  NodeValue_rescaled[measnode[im]], NodeError_rescaled[measnode[im]],
		  chisquare_meas_rescaled[im], chisquare_node, chisquare_tot_rescaled,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g \n\n", 
	       chisquare_tot_rescaled, nchisquare_tot_rescaled,
	       nbase_u, nchisquare_tot_rescaled - nbase_u, TMath::Prob(chisquare_tot_rescaled, nchisquare_tot_rescaled - nbase_u));
  //
  cout << Form("%s = (%10.6f +- %10.6f)%% (%10.6f sigma); %s = (%10.6f +- %10.6f)%% (%10.6f sigma)\n\n",
	       nodetitle[N_GAMMAALL].data(),
	       NodeValue_rescaled[N_GAMMAALL]*100.,
	       uconstrain?0:NodeError_rescaled[N_GAMMAALL]*100.,
	       uconstrain?0:(1-NodeValue_rescaled[N_GAMMAALL])/NodeError_rescaled[N_GAMMAALL],
	       nodetitle[N_GAMMA110].data(),
	       NodeValue_rescaled[N_GAMMA110]*100.,
	       NodeError_rescaled[N_GAMMA110]*100.,
	       NodeValue_rescaled[N_GAMMA110]/NodeError_rescaled[N_GAMMA110]);
  //
  cout << "Summary of fit with [errors rescaled in a Ad-Hoc style for kkk only]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE             SEED      FITVAL         FITERR     SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %6.3f %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),baseseed[ibase],
		 basevalue_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile_rescaled[2];
  for (int p=0;p<2;++p){
    if (p==0) measfile_rescaled[0]=fopen(Form("combos_measurements_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile_rescaled[1]=fopen(Form("alucomb_measurements_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile_rescaled[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_none_rescaled, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror_rescaled, corrmat,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile_rescaled[2];
  for (int p=0;p<2;++p){
    if (p==0) avefile_rescaled[0]=fopen(Form("combos_average_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile_rescaled[1]=fopen(Form("alucomb_average_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile_rescaled[p], "INCLUDE combos_measurements_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile_rescaled[p], "INCLUDE alucomb_measurements_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile_rescaled[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_none_rescaled, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, node_is_base,
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  //  Prepare to obtain Fit Values for non-weak measurements [rescaled]
  //
  double basevalue_fit_noweak_rescaled[nbase];
  double baseerror_fit_noweak_rescaled[nbase];
  double** basecov_fit_noweak_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_noweak_rescaled[i] = new double[nbase];
  double** basecorr_fit_noweak_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_noweak_rescaled[i] = new double[nbase];
  double chisquared_noweak_rescaled;
  //
  // Define weak nodes
  //
  double nsig_rescaled[200];
  int weak_rescaled[200];
  for (i=0;i<nmeas;++i){
    if (nchisquare_meas_rescaled[i] > 0 ) {
      nsig_rescaled[i] = measerror[i]/(TMath::Sqrt(nchisquare_meas_rescaled[i])*NodeError_rescaled[measnode[i]]);
    } else {
      nsig_rescaled[i] = 999999;
    }
    weak_rescaled[i] = (icorrj[i].size()==1) && (nsig_rescaled[i] > 3.); // only for uncorrelated measurements
  }
  //
  // COMBINE
  //
  combine(uconstrain,
	  nmeas, weak_rescaled, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix_rescaled,
	  node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_noweak_rescaled, // output
	  baseerror_fit_noweak_rescaled, // output
	  basecov_fit_noweak_rescaled,   // output
	  basecorr_fit_noweak_rescaled,  // output
	  chisquared_noweak_rescaled     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_rescaled,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_noweak_rescaled_temp = chisquared_noweak_rescaled;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_noweak_rescaled[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE
    //
    combine(uconstrain,
	    nmeas, weak_rescaled, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix_rescaled,
	    node_num_npar, node_den_npar, node_quan, node_parm, node_is_base,
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_noweak_rescaled, // output
	    baseerror_fit_noweak_rescaled, // output
	    basecov_fit_noweak_rescaled,   // output
	    basecorr_fit_noweak_rescaled,  // output
	    chisquared_noweak_rescaled     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared_noweak_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_rescaled,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_noweak_rescaled - chisquared_noweak_rescaled_temp) < 1e-3) break;
    chisquared_noweak_rescaled_temp = chisquared_noweak_rescaled;
  }
  //
  cout << endl << "Comparison of Results from fit with [non-weak measurements only] and [errors rescaled in a Ad-Hoc style for kkk only] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE      ORIG_FITVAL    ORIG_FITERR    SCAL_FITVAL    SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %14.7f %6.3f  %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_noweak_rescaled[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_noweak_rescaled[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_noweak_rescaled[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_noweak_rescaled[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_noweak_rescaled[nnode];     // value of each node
  vector<double> NodeValue_part_noweak_rescaled[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_noweak_rescaled, // output
		   NodeValue_den_noweak_rescaled, // output
		   NodeValue_noweak_rescaled,     // output
		   NodeValue_part_noweak_rescaled // output
		   );
  //
  double NodeError_noweak_rescaled[nnode];
  double** NodeErrorMatrix_noweak_rescaled = new double*[nnode]; for (i=0;i<nnode;++i) NodeErrorMatrix_noweak_rescaled[i]= new double[nnode];
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_noweak_rescaled[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_noweak_rescaled[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_noweak_rescaled[jnode].at(jpar);
	  NodeErrorMatrix_noweak_rescaled[inode][jnode] += (iquan>nbase || jquan>nbase) ? 0 :
	    ipartial*jpartial*baseerror_fit_noweak_rescaled[iquan-1]*baseerror_fit_noweak_rescaled[jquan-1]*basecorr_fit_noweak_rescaled[iquan-1][jquan-1];
	}
      }
    }
    if (NodeErrorMatrix_noweak_rescaled[inode][inode]<0) NodeErrorMatrix_noweak_rescaled[inode][inode] = 0;
    //
    NodeError_noweak_rescaled[inode] = TMath::Sqrt(TMath::Max(0.,NodeErrorMatrix_noweak_rescaled[inode][inode]));
  }
  //
  int   nchisquare_tot_noweak_rescaled=0; 
  double chisquare_tot_noweak_rescaled=0; 
  double chisquare_meas_noweak_rescaled[200]; // chi2 per mesurement
  int   nchisquare_meas_noweak_rescaled[200];  for (im=0;im<nmeas;++im) nchisquare_meas_noweak_rescaled[im]=0;// nodes per mesurement     
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_noweak_rescaled[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_rescaled[jm]==weakcompare) continue;
	  chisquare_meas_noweak_rescaled[im]+= 
	    (measvalue[im]-NodeValue_noweak_rescaled[measnode[im]])*InvMeasErrorMatrix_rescaled[im][jm]*(measvalue[jm]-NodeValue_noweak_rescaled[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_noweak_rescaled[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_noweak_rescaled+=1;
	chisquare_tot_noweak_rescaled+=chisquare_node;
	nchisquare_meas_noweak_rescaled[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror_rescaled[im], 
		  NodeValue_noweak_rescaled[measnode[im]], NodeError_noweak_rescaled[measnode[im]],
		  chisquare_meas_noweak_rescaled[im], chisquare_node, chisquare_tot_noweak_rescaled,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g \n\n", 
	       chisquare_tot_noweak_rescaled, nchisquare_tot_noweak_rescaled,
	       nbase_u, nchisquare_tot_noweak_rescaled - nbase_u, TMath::Prob(chisquare_tot_noweak_rescaled, nchisquare_tot_noweak_rescaled - nbase_u));
  //
  cout << Form("%s = (%10.6f +- %10.6f)%% (%10.6f sigma); %s = (%10.6f +- %10.6f)%% (%10.6f sigma)\n\n",
	       nodetitle[N_GAMMAALL].data(),
	       NodeValue_noweak_rescaled[N_GAMMAALL]*100.,
	       uconstrain?0:NodeError_noweak_rescaled[N_GAMMAALL]*100.,
	       uconstrain?0:(1-NodeValue_noweak_rescaled[N_GAMMAALL])/NodeError_noweak_rescaled[N_GAMMAALL],
	       nodetitle[N_GAMMA110].data(),
	       NodeValue_noweak_rescaled[N_GAMMA110]*100.,
	       NodeError_noweak_rescaled[N_GAMMA110]*100.,
	       NodeValue_noweak_rescaled[N_GAMMA110]/NodeError_noweak_rescaled[N_GAMMA110]);
  //
  cout << "Summary of fit with [non-weak measurements only] and [errors rescaled in a Ad-Hoc style for kkk only]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE             SEED      FITVAL         FITERR     SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %14.7f %14.7f %14.7f %6.3f %-60s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenodename[ibase].data(),baseseed[ibase],
		 basevalue_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile_noweak_rescaled[2];
  for (int p=2;p<2;++p){
    if (p==0) measfile_noweak_rescaled[0]=fopen(Form("combos_measurements_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile_noweak_rescaled[1]=fopen(Form("alucomb_measurements_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile_noweak_rescaled[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_rescaled, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror_rescaled, corrmat,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile_noweak_rescaled[2];
  for (int p=2;p<2;++p){
    if (p==0) avefile_noweak_rescaled[0]=fopen(Form("combos_average_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile_noweak_rescaled[1]=fopen(Form("alucomb_average_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile_noweak_rescaled[p], "INCLUDE combos_measurements_noweak_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile_noweak_rescaled[p], "INCLUDE alucomb_measurements_noweak_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile_noweak_rescaled[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_rescaled, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, node_is_base,
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  // Print final values used from fit with [errors rescaled in a Ad-Hoc style for kkk only]
  //
  cout << endl << "Final Results from fit with [errors rescaled in a Ad-Hoc style for kkk only] w.r.t original fit:" << endl;
  cout << Form("\\begin{tabular}{l c}\n");
  for (i=1;i<=nbase;++i) {
    for (ibase=0;ibase<nbase;++ibase) {
      if (i!=baseorder[ibase]) continue;
#if defined USING_NBASE39
      if (ibase==M_GAMMA128) cout << Form("%-70s &                     \\\\ \n","$K^- \\phi \\nu_\\tau (\\phi \\to KK)$");
#endif
      if (basevalue_fit_rescaled[ibase]*100>10) {
	cout << Form("%-70s & %6.3f $\\pm$ %6.3f \\\\ \n",baselatex[ibase].data(),basevalue_fit_rescaled[ibase]*100,baseerror_fit_rescaled[ibase]*100);
      } else {
	cout << Form("%-70s & %5.3f $\\pm$ %5.3f \\\\ \n",baselatex[ibase].data(),basevalue_fit_rescaled[ibase]*100,baseerror_fit_rescaled[ibase]*100);
      }
    }
  }
  cout << Form("\\end{tabular}\n");
  cout << endl;
  for (inode=0;inode<nnode;++inode) {
    cout << Form("inode=%4d %-8s %-12s %-108s Fit = %8.4e +- %8.4e\n",
		 inode, a_nodename[inode], nodegammaname[inode].data(), nodetitle[inode].data(), NodeValue_rescaled[inode], NodeError_rescaled[inode]);
  }
  cout << endl;
  //
  // Set Global Parameters
  //
  GlobalPar[P_m_e] =  m_e; e_GlobalPar[P_m_e] = e_m_e;
  GlobalPar[P_m_mu] = m_mu; e_GlobalPar[P_m_mu] = e_m_mu;
  GlobalPar[P_tau_mu] = tau_mu; e_GlobalPar[P_tau_mu] = e_tau_mu;
  GlobalPar[P_m_tau] = m_tau; e_GlobalPar[P_m_tau] =  e_m_tau;
  GlobalPar[P_tau_tau] = tau_tau; e_GlobalPar[P_tau_tau] = e_tau_tau;
  GlobalPar[P_m_pim] = m_pim; e_GlobalPar[P_m_pim] = e_m_pim;
  GlobalPar[P_tau_pim] = tau_pim; e_GlobalPar[P_tau_pim] = e_tau_pim;
  GlobalPar[P_BR_PimToMumNu] = BR_PimToMumNu; e_GlobalPar[P_BR_PimToMumNu] = e_BR_PimToMumNu;
  GlobalPar[P_Delta_TauToPim_over_PimToMu] = Delta_TauToPim_over_PimToMu; e_GlobalPar[P_Delta_TauToPim_over_PimToMu] = e_Delta_TauToPim_over_PimToMu;
  GlobalPar[P_m_km] = m_km; e_GlobalPar[P_m_km] = e_m_km;
  GlobalPar[P_tau_km] = tau_km; e_GlobalPar[P_tau_km] = e_tau_km;
  GlobalPar[P_BR_KmToMumNu] = BR_KmToMumNu; e_GlobalPar[P_BR_KmToMumNu] = e_BR_KmToMumNu;
  GlobalPar[P_Delta_TauToKm_over_KmToMu] = Delta_TauToKm_over_KmToMu; e_GlobalPar[P_Delta_TauToKm_over_KmToMu] = e_Delta_TauToKm_over_KmToMu;
  GlobalPar[P_m_W] = m_W; e_GlobalPar[P_m_W] = e_m_W;
  GlobalPar[P_alpha_val] = alpha_val; e_GlobalPar[P_alpha_val] = alpha_err;
  GlobalPar[P_Delta_mu_gamma] = Delta_mu_gamma; e_GlobalPar[P_Delta_mu_gamma] = e_Delta_mu_gamma;
  GlobalPar[P_Delta_tau_gamma] = Delta_tau_gamma; e_GlobalPar[P_Delta_tau_gamma] = e_Delta_tau_gamma;
  GlobalPar[P_GMu] = GMu; e_GlobalPar[P_GMu] = e_GMu;
  GlobalPar[P_hbar] = hbar; e_GlobalPar[P_hbar] = e_hbar;
  GlobalPar[P_SEW] = SEW; e_GlobalPar[P_SEW] = e_SEW;
  GlobalPar[P_fpi] = fpi; e_GlobalPar[P_fpi] = e_fpi;
  GlobalPar[P_fK] = fK; e_GlobalPar[P_fK] = e_fK;
  GlobalPar[P_fKfpi] = fKfpi; e_GlobalPar[P_fKfpi] = e_fKfpi;
  GlobalPar[P_Vud] = Vud; e_GlobalPar[P_Vud] = e_Vud;
  GlobalPar[P_Delta_Kpi] = Delta_Kpi; e_GlobalPar[P_Delta_Kpi] = e_Delta_Kpi;
  GlobalPar[P_Delta_Rth] = Delta_Rth; e_GlobalPar[P_Delta_Rth] = e_Delta_Rth;
  //
  // Print Results 
  //
  FILE *resfile=fopen(Form("readpdg_%s%s.results",sconstrain.data(),salephhcorr.data()),"w");
  for (i=0;i<3;++i) {
    if (i==0) {
      fprintf (resfile, "Results from original fit :\n\n");
      calc_results(resfile,
		   basetitle, basevalue_fit, baseerror_fit, basecov_fit, basecorr_fit, 
		   NodeValue, NodeError, NodeErrorMatrix);
      fprintf (resfile, "\n==============================================\n");
    } else if (i==1) {
      fprintf (resfile, "Results from fit [errors inflated with PDG-style scale factors] :\n\n");
      calc_results(resfile,
		   basetitle, basevalue_fit, baseerror_fit_noweak_scaled, basecov_fit_noweak_scaled, basecorr_fit_noweak_scaled, 
		   NodeValue, NodeError_noweak_scaled, NodeErrorMatrix_noweak_scaled);
      fprintf (resfile, "\n==============================================\n");
    } else if (i==2) {
      fprintf (resfile, "Results from fit [errors rescaled in a Ad-Hoc style for kkk only] :\n\n");
      calc_results(resfile,
		   basetitle, basevalue_fit_rescaled, baseerror_fit_rescaled, basecov_fit_rescaled, basecorr_fit_rescaled, 
		   NodeValue_rescaled, NodeError_rescaled, NodeErrorMatrix_rescaled);
      fprintf (resfile, "\n==============================================\n");
    }
  }
  fclose(resfile);
  // Clean Up
  //
  delete [] NodeErrorMatrix_noweak_rescaled;
  delete [] basecorr_fit_noweak_rescaled;
  delete [] basecov_fit_noweak_rescaled;
  delete [] NodeErrorMatrix_rescaled;
  delete [] basecorr_fit_rescaled;
  delete [] basecov_fit_rescaled;
  delete [] NodeErrorMatrix_noweak_scaled;
  delete [] basecorr_fit_noweak_scaled;
  delete [] basecov_fit_noweak_scaled;
  delete [] corrmat_scaled_noweak;
  delete [] corrmat_scaled;
  delete [] NodeErrorMatrix_noweak;
  delete [] basecorr_fit_noweak;
  delete [] basecov_fit_noweak;
  delete [] NodeErrorMatrix;
  delete [] basecorr_fit;
  delete [] basecov_fit;
  delete [] corrmat;
  delete [] a_nodename;
}
