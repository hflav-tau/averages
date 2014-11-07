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
#include <cstdlib>
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
// ----------------------------------------------------------------------
int main(int argc, char* argv[]){
  int ifile;
  char firstch;
  string dummy;
  //
  // INPUTS:
  //
  string fname_common = "../Common/Parameters.input" ;
  //
  vector<string> fnames;
  fnames.push_back("../TauFit/babar/gamma3by5.input");
  fnames.push_back("../TauFit/babar/gamma9by5.input");
  fnames.push_back("../TauFit/babar/gamma10by5.input");
  fnames.push_back("../TauFit/babar/gamma16.input");
  fnames.push_back("../TauFit/babar/gamma35.input");
  fnames.push_back("../TauFit/babar/gamma40.input");
  fnames.push_back("../TauFit/babar/gamma60.input");
  fnames.push_back("../TauFit/babar/gamma85.input");
  fnames.push_back("../TauFit/babar/gamma93.input");
  fnames.push_back("../TauFit/babar/gamma103.input");
  fnames.push_back("../TauFit/babar/gamma96.input");
  fnames.push_back("../TauFit/babar/gamma136.input");
  fnames.push_back("../TauFit/babar/gamma128.input");
  fnames.push_back("../TauFit/belle/gamma13.input");
  fnames.push_back("../TauFit/belle/gamma35.input");
  fnames.push_back("../TauFit/belle/gamma60.input");
  fnames.push_back("../TauFit/belle/gamma85.input");
  fnames.push_back("../TauFit/belle/gamma93.input");
  fnames.push_back("../TauFit/belle/gamma126.input");
  fnames.push_back("../TauFit/belle/gamma128.input");
  fnames.push_back("../TauFit/belle/gamma130.input");
  fnames.push_back("../TauFit/belle/gamma132.input");
  fnames.push_back("../TauFit/belle/gamma96.input");
  //  fnames.push_back("../TauFit/belle/gamma115.input");
  //
  int ipar;
  string par_name[100];
  double par_val[100], par_poserr[100], par_negerr[100];
  //
  ifstream ifs(fname_common.data());
  if (!ifs.good()) {
    cout << Form("Cannot open input file : %s\n", fname_common.data()); 
    exit(1);
  }
  //
  ipar=-1;
  int nlines=-1;
  while(ifs.good()) {
    ++nlines;
    firstch=' '; ifs.get(firstch) ;  ifs.putback(firstch) ;
    if (firstch=='*'||firstch=='\n') { // Skip such lines
      ifs.ignore(256,'\n') ;
    } else {
      ifs >> dummy;
      if (strcmp(dummy.data(),"PARAMETERS")==0) {
	ifs.ignore(256,'\n') ; 
      } else {
	++ipar;
	par_name[ipar].assign(dummy);
	ifs >> par_val[ipar] >> par_poserr[ipar] >> par_negerr[ipar];
	ifs.ignore(256,'\n') ;
      }
    }
  }
  ifs.close();
  int npar=ipar;
  //
  cout << Form("Read %d lines, %d parameters from : %s\n", nlines,npar,fname_common.data());
  for (ipar=0;ipar<npar;++ipar){
    cout << "ipar = " << ipar << " " << par_name[ipar] << " " << par_val[ipar] << " " << par_poserr[ipar] << " " << par_negerr[ipar] << endl;
  }
  //
  //
  int imeas=-1;
  int jmeas, kpar, kpar_sv, esys, esys_sv;
  vector<string> nsys_name;
  string meas_name[100], meas_exp[100], meas_gamma[100], meas_pub[100];
  string meas_corr_exp[100][100], meas_corr_gamma[100][100], meas_corr_pub[100][100], meas_par_name[100][100], meas_esys_name[100][100];
  double meas_val[100], meas_stat[100], meas_syst[100], meas_corr_val[100][100], meas_val_adj[100], meas_err_adj[100];
  double meas_par_value[100][100], meas_par_poserr[100][100], meas_par_negerr[100][100], meas_esys_val[100][100], meas_esys_val_adj[100][100];
  double meas_cov[100][100], meas_corr[100][100];
  int meas_corr_num[100], meas_corrij[100][100], meas_par_num[100], meas_sys_num[100];
  for (ifile=0;ifile<fnames.size(); ++ifile) {
    vector<string> vdummy;
    ifs.open(fnames[ifile].data(),ifstream::in);
    if (!ifs.good()) {
      cout << Form("Cannot open input file : %s\n", fnames[ifile].data()); 
      exit(1);
    } else {
      while (ifs.good()){
	firstch=' '; ifs.get(firstch); ifs.putback(firstch) ;
	if (firstch=='*' || firstch=='\n') {// Skip such lines
	  ifs.ignore(256,'\n'); 
	  continue;
	}
	ifs >> dummy;
	if (strcmp(dummy.data(),"*")==0||strcmp(dummy.data(),"!")==0) { // Skip rest of the line
	  ifs.ignore(256,'\n'); 
	  continue;
	}
	if (
	    (strcmp(dummy.data(),"statistical")==0) ||
	    (strcmp(dummy.data(),"systematic")==0) ||
	    (strcmp(dummy.data(),"stat")==0) ||
	    (strcmp(dummy.data(),"syst")==0)
	    ) {/*skip*/} else { // Keep the rest
	  vdummy.push_back(dummy);
	}
      }
    }
    ifs.close();
    //    cout << Form("Read from file: %s\n",fnames[ifile].data());
    //
    bool kpar_start=false;
    bool esys_start=false;
    for (vector<string>::iterator vsi=vdummy.begin();vsi!=vdummy.end();++vsi) {
      dummy=*vsi;
      int keyword = 0;
      if (
	  (strcmp(dummy.data(),"BEGIN")==0) ||
	  (strcmp(dummy.data(),"PARAMETERS")==0) ||
	  (strcmp(dummy.data(),"MEASUREMENT")==0) ||
	  (strcmp(dummy.data(),"STAT_CORR_WITH")==0) ||
	  (strcmp(dummy.data(),"DATA")==0) ||
	  (strcmp(dummy.data(),"END")==0)
	  ) {
	keyword = 1;
      }
      if (strcmp(dummy.data(),"BEGIN")==0) {
	++imeas;
	jmeas=-1;
	kpar=-2;
	esys=-2;
	meas_exp[imeas]=*(vsi+1);
	meas_gamma[imeas]=*(vsi+2); meas_gamma[imeas].erase(0,5);
	meas_pub[imeas]=*(vsi+3);
      }
      if (strcmp(dummy.data(),"MEASUREMENT")==0){
	meas_name[imeas]=*(vsi+1);
	meas_val[imeas]=strtod((*(vsi+4)).data(),0);
	meas_stat[imeas]=strtod((*(vsi+5)).data(),0);
	meas_syst[imeas]=strtod((*(vsi+6)).data(),0);
      }
      if (strcmp(dummy.data(),"STAT_CORR_WITH")==0){
	++jmeas;
	meas_corr_exp[imeas][jmeas]=*(vsi+1);
	meas_corr_gamma[imeas][jmeas]=*(vsi+2); meas_corr_gamma[imeas][jmeas].erase(0,5);
	meas_corr_pub[imeas][jmeas]=*(vsi+3);
	meas_corr_val[imeas][jmeas]=strtod((*(vsi+4)).data(),0);
      }
      if (strcmp(dummy.data(),"PARAMETERS")==0){
	kpar_start=true;
      } else {
	if (kpar_start && keyword) kpar_start=false; 
      }
      if (kpar_start) ++kpar;
      if (kpar>-1 && kpar!=kpar_sv) {
	if (kpar%4==0) meas_par_name[imeas][kpar/4] = *vsi;
	if (kpar%4==1) meas_par_value[imeas][kpar/4] =  strtod((*vsi).data(),0);
	if (kpar%4==2) meas_par_poserr[imeas][kpar/4] =  strtod((*vsi).data(),0);
	if (kpar%4==3) meas_par_negerr[imeas][kpar/4] =  strtod((*vsi).data(),0);
      }
      if (strcmp(dummy.data(),"DATA")==0 && strcmp((*(vsi+1)).substr(0,2).data(),"m_")!=0) {
	//	cout << *(vsi+1) << " " << strcmp((*(vsi+1)).substr(0,2).data(),"m_") << endl;
	esys_start=true;
      } else {
	if (esys_start && keyword) esys_start=false;
      }
      if (esys_start) ++esys;
      //
      if (esys>-1 && esys!=esys_sv) {
	if (esys%2==0) {
	  meas_esys_name[imeas][esys/2] = *vsi; 
	  string stemp=*vsi;
	  size_t found=stemp.find("%");
	  if (found!=string::npos) stemp.erase(found);
	  nsys_name.push_back(stemp);
	}
	if (esys%2==1) meas_esys_val_adj[imeas][esys/2] = meas_esys_val[imeas][esys/2] = strtod((*vsi).data(),0); 
      }
      //
      kpar_sv = kpar;
      esys_sv = esys;
      /*
      cout  << " keyword = " << keyword << " dummy = " << dummy // << " substr  = " << dummy.substr(0,2) 
	    << " kpar = " << kpar << " kpar_sv = " << kpar_sv 
	    << " esys = " << esys << " esys_sv = " << esys_sv 
	    << endl;
      */
    }
    //
    meas_corr_num[imeas] = jmeas+1;
    meas_par_num[imeas] = (kpar+1)/4;
    meas_sys_num[imeas] = (esys+1)/2;
    //
  }
  //
  sort(nsys_name.begin(),nsys_name.end());
  vector<string>::iterator new_end=unique(nsys_name.begin(),nsys_name.end());
  int nsys=nsys_name.size();
  cout << "nsys = " << nsys << " :: " ; for (esys=0;esys<nsys;++esys) {cout << esys << " " << nsys_name[esys] << " , "; } cout << endl << endl;
  nsys=new_end-nsys_name.begin();
  cout << "nsys = " << nsys << " :: " ; for (esys=0;esys<nsys;++esys) {cout << esys << " " << nsys_name[esys] << " , "; } cout << endl << endl;
  //
  int nmeas=imeas+1;
  for (imeas=0;imeas<nmeas;++imeas){
    cout << endl;
    cout << Form("Read from file: %s\n",fnames[imeas].data());
    cout << "imeas = " << imeas 
	 << " name = " << meas_name[imeas] 
	 << " exp = " <<  meas_exp[imeas] 
	 << " gamma = " << meas_gamma[imeas] 
	 << " pub = " << meas_pub[imeas]  
	 << " val = " << meas_val[imeas] 
	 << " stat = " << meas_stat[imeas] 
	 << " syst = " << meas_syst[imeas] 
	 << " corr_n = " << meas_corr_num[imeas]
	 << " par_n = " << meas_par_num[imeas] 
	 << " sys_n = " << meas_sys_num[imeas]
	 << endl;
    //
    for (esys=0;esys<meas_sys_num[imeas];++esys) {
      size_t found=meas_esys_name[imeas][esys].find("%");
      if (found!=string::npos) {
	meas_esys_name[imeas][esys].erase(found);
	meas_esys_val[imeas][esys]*=meas_val[imeas]/100.;	
	meas_esys_val_adj[imeas][esys]*=meas_val[imeas]/100.;	
      }
      cout << "esys = " << esys << " name = " << meas_esys_name[imeas][esys] << " val = " << meas_esys_val[imeas][esys] << " adj = " << meas_esys_val_adj[imeas][esys] << endl;
    }
    //
    meas_val_adj[imeas] = meas_val[imeas];
    for (kpar=0;kpar<meas_par_num[imeas];++kpar) {
      int this_esys=-1;
      for (esys=0;esys<meas_sys_num[imeas];++esys) {
	if (strcmp(meas_par_name[imeas][kpar].data(),meas_esys_name[imeas][esys].data())==0) {
	  this_esys=esys;
	  break;
	}
      }
      if (this_esys==-1) continue;
      for (ipar=0;ipar<npar;++ipar){
	if (strcmp(meas_par_name[imeas][kpar].data(),par_name[ipar].data())==0) {
	  meas_val_adj[imeas]+= ((meas_par_value[imeas][kpar] - par_val[ipar])/(0.5*(meas_par_poserr[imeas][kpar] - meas_par_negerr[imeas][kpar]))) * fabs(meas_esys_val[imeas][this_esys]);
	  meas_esys_val_adj[imeas][this_esys] = meas_esys_val[imeas][this_esys] * (0.5*(par_poserr[ipar] - par_negerr[ipar])) / (0.5*(meas_par_poserr[imeas][kpar] - meas_par_negerr[imeas][kpar])) ;
	}
      }
      cout << "kpar = " << kpar << " name = " << meas_par_name[imeas][kpar] << " val = " << meas_par_value[imeas][kpar] << " poserr = " << meas_par_poserr[imeas][kpar] << " negerr = " << meas_par_negerr[imeas][kpar] << " meas_val_adj = " << meas_val_adj[imeas] << " esys = " << this_esys << " esys_val_adj = " << meas_esys_val_adj[imeas][this_esys] << endl;
    }
    meas_err_adj[imeas] = pow(meas_stat[imeas],2);
    for (esys=0;esys<meas_sys_num[imeas];++esys) {
      meas_err_adj[imeas]+= pow(meas_esys_val_adj[imeas][esys],2);
    }
    meas_err_adj[imeas] = sqrt(max(0.,meas_err_adj[imeas]));
    cout << Form("adjusted value = %20.10fg +- %20.10f\n",meas_val_adj[imeas],meas_err_adj[imeas]) ;
    //
    for (jmeas=0;jmeas<meas_corr_num[imeas];++jmeas) {
      for (int ijmeas=0;ijmeas<nmeas;++ijmeas) {
	if ((strcmp(meas_exp[ijmeas].data(),meas_corr_exp[imeas][jmeas].data())==0) &&
	    (strcmp(meas_gamma[ijmeas].data(),meas_corr_gamma[imeas][jmeas].data())==0) &&
	    (strcmp(meas_pub[ijmeas].data(),meas_corr_pub[imeas][jmeas].data())==0)) {
	  meas_corrij[imeas][jmeas]=ijmeas;
	  break;
	}
      }
      cout << "jmeas = " << jmeas << " exp = " << meas_corr_exp[imeas][jmeas] << " gamma = " << meas_corr_gamma[imeas][jmeas] << " pub = " << meas_corr_pub[imeas][jmeas] << " corrij = " << meas_corrij[imeas][jmeas] << " val = " << meas_corr_val[imeas][jmeas] << endl;
    }
    //
    for (jmeas=0;jmeas<imeas;++jmeas) {
      meas_cov[imeas][jmeas]=0;
      for (int ijmeas=0;ijmeas<meas_corr_num[imeas];++ijmeas) {
	if (jmeas==meas_corrij[imeas][ijmeas]) {
	  meas_cov[imeas][jmeas]+=meas_corr_val[imeas][ijmeas]*meas_stat[imeas]*meas_stat[jmeas];
	  cout << "imeas = " << imeas << " stat = " << meas_stat[imeas] << " jmeas = " << jmeas << " stat = " << meas_stat[jmeas] << " ijmeas = " << ijmeas << " corr = " << meas_corr_val[imeas][ijmeas] << " cov = " << meas_cov[imeas][jmeas] << endl;
	}
      }
      for (esys=0;esys<nsys;++esys){	  
	int i_esys=-1;
	for (i_esys=0;i_esys<meas_sys_num[imeas];++i_esys) {
	  if (strcmp(meas_esys_name[imeas][i_esys].data(),nsys_name[esys].data())==0) {
	    //	    cout << meas_esys_name[imeas][i_esys] << endl;
	    break;
	  }
	}
	int j_esys=-1;
	for (j_esys=0;j_esys<meas_sys_num[jmeas];++j_esys) {
	  if (strcmp(meas_esys_name[jmeas][j_esys].data(),nsys_name[esys].data())==0) {
	    //	    cout << meas_esys_name[jmeas][j_esys] << endl;
	    break;
	  }
	}
	if (i_esys<meas_sys_num[imeas] && j_esys<meas_sys_num[jmeas]){
	  cout << " imeas = " << imeas << " name = " << meas_name[imeas] << " i_esys = " << i_esys << " name = " << meas_esys_name[imeas][i_esys]
	       << " jmeas = " << jmeas << " name = " << meas_name[jmeas] << " j_esys = " << j_esys << " name = " << meas_esys_name[jmeas][j_esys] << endl;
	  meas_cov[imeas][jmeas]+=meas_esys_val_adj[imeas][i_esys] * meas_esys_val_adj[jmeas][j_esys];
	}
      }
      meas_corr[imeas][jmeas] = meas_cov[imeas][jmeas] / (meas_err_adj[imeas]*meas_err_adj[jmeas]) ;
      if ( meas_cov[imeas][jmeas] != 0) {
	cout << "meas_cov["<<imeas<<"]["<<jmeas<<"] = " << Form("%20.10f",meas_cov[imeas][jmeas]) << " " 
	     << "meas_cor["<<imeas<<"]["<<jmeas<<"] = " << Form("%20.10f",meas_corr[imeas][jmeas]) << " "
	     << endl;
      }
    }
    //
  }
  for (imeas=0;imeas<nmeas;++imeas){ 
    for (jmeas=0;jmeas<nmeas;++jmeas) {
      if (imeas==jmeas){
	meas_cov[imeas][jmeas] = pow(meas_err_adj[imeas],2);
	meas_corr[imeas][jmeas]=1;
      }
      if (jmeas>imeas) {
	meas_cov[imeas][jmeas] = meas_cov[jmeas][imeas];
	meas_corr[imeas][jmeas]= meas_corr[jmeas][imeas];
      }
    }
  }
  //
  FILE *thisfile = fopen("readBB.out","w");
  fprintf (thisfile,"\nbabar belle measurements \n\n");
  for (imeas=0;imeas<nmeas;++imeas){ 
    string istr= meas_exp[imeas]+".Gamma"+meas_gamma[imeas]+"."+meas_pub[imeas];
    fprintf (thisfile,"    %-40s\n", istr.data());
    fprintf (thisfile,"val       %20.10f\n", meas_val_adj[imeas]);
    fprintf (thisfile,"err       %20.10f\n", meas_err_adj[imeas]);
  }
  fprintf (thisfile,"\nbabar belle correlation\n\n");
  for (imeas=0;imeas<nmeas;++imeas){ 
    string istr= meas_exp[imeas]+".Gamma"+meas_gamma[imeas]+"."+meas_pub[imeas];
    fprintf (thisfile,"                           %-40s\n", istr.data());
    for (jmeas=0;jmeas<nmeas;++jmeas){ 
      string jstr= meas_exp[jmeas]+".Gamma"+meas_gamma[jmeas]+"."+meas_pub[jmeas];
      fprintf (thisfile,"%-30s  %20.10f\n", jstr.data(), meas_corr[imeas][jmeas]);
    }
  }
  fclose(thisfile);
  //
  //less all_node_def.txt | awk '{if ($1=="*") print $0}' | sort -k2 -g | sort -u | sed -e 's/* S035/S035/g' | awk -F" : " '{print "++inode; node_name[inode] = \""$1"\"; node_title[inode] = \""$2"\"; node_gamma[inode] = \""$3"\";"}'
  //
  int inode=-1;
  string node_name[200], node_title[200], node_gamma[200];
  ++inode; node_name[inode] = "S035B20"; node_title[inode] = "G(K- eta nu(tau)) / G(total)"; node_gamma[inode] = "Gamma128";
  ++inode; node_name[inode] = "S035B21"; node_title[inode] = "G(h- 2pi0 nu(tau) (ex. K0)) / G(h- pi0 nu(tau))"; node_gamma[inode] = "Gamma19by13";
  ++inode; node_name[inode] = "S035B22"; node_title[inode] = "G(h- 3pi0 nu(tau)) / G(h- pi0 nu(tau))"; node_gamma[inode] = "Gamma26by13";
  ++inode; node_name[inode] = "S035B23"; node_title[inode] = "G(h- 4pi0 nu(tau) (ex. K0, eta)) / G(total)"; node_gamma[inode] = "Gamma30";
  ++inode; node_name[inode] = "S035B25"; node_title[inode] = "G(h- h- h+ 2pi0 nu(tau) (ex. K0)) / G(h- h- h+ >=0 neutrals >=0 K(L)0 nu(tau))"; node_gamma[inode] = "Gamma76by54";
  ++inode; node_name[inode] = "S035B26"; node_title[inode] = "G(h- omega pi0 nu(tau)) / G(h- h- h+ >=0 neutrals >=0 K(L)0 nu(tau))"; node_gamma[inode] = "Gamma152by54";
  ++inode; node_name[inode] = "S035B27"; node_title[inode] = "G(h- omega pi0 nu(tau)) / G(h- h- h+ 2pi0 nu(tau) (ex. K0))"; node_gamma[inode] = "Gamma152by76";
  ++inode; node_name[inode] = "S035B29"; node_title[inode] = "G(K- pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma16";
  ++inode; node_name[inode] = "S035B30"; node_title[inode] = "G(K- 2pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma23";
  ++inode; node_name[inode] = "S035B31"; node_title[inode] = "G(K- 3pi0 nu(tau) (ex. K0, eta)) / G(total)"; node_gamma[inode] = "Gamma28";
  ++inode; node_name[inode] = "S035B32"; node_title[inode] = "G(Kbar0 pi- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma35";
  ++inode; node_name[inode] = "S035B33"; node_title[inode] = "G(Kbar0 pi- pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma40";
  ++inode; node_name[inode] = "S035B34"; node_title[inode] = "G(K- pi0 K0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma42";
  ++inode; node_name[inode] = "S035B37"; node_title[inode] = "G(pi- K- K+ >=0 neutrals nu(tau)) / G(total)"; node_gamma[inode] = "Gamma92";
  ++inode; node_name[inode] = "S035B43"; node_title[inode] = "G(K(S)0 (particles)- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma33";
  ++inode; node_name[inode] = "S035B45"; node_title[inode] = "G((5pi)- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma106";
  ++inode; node_name[inode] = "S035B51"; node_title[inode] = "G(pi- K0 Kbar0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma46";
  ++inode; node_name[inode] = "S035B53"; node_title[inode] = "G(h- h- h+ pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma66";
  ++inode; node_name[inode] = "S035B54"; node_title[inode] = "G(h- h- h+ pi0 nu(tau) (ex. K0, omega)) / G(total)"; node_gamma[inode] = "Gamma67";
  ++inode; node_name[inode] = "S035B55"; node_title[inode] = "G(pi- 2pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma20";
  ++inode; node_name[inode] = "S035B56"; node_title[inode] = "G(pi- 3pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma27";
  ++inode; node_name[inode] = "S035B57"; node_title[inode] = "G(h- h- h+ 3pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma78";
  ++inode; node_name[inode] = "S035B58"; node_title[inode] = "G(h- pi0 omega nu(tau)) / G(total)"; node_gamma[inode] = "Gamma152";
  ++inode; node_name[inode] = "S035B59"; node_title[inode] = "G(h- h- h+ 2pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma76";
  ++inode; node_name[inode] = "S035B60"; node_title[inode] = "G(K- Kstar(892)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma115";
  ++inode; node_name[inode] = "S035B62"; node_title[inode] = "G(h- h- h+ nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma57";
  ++inode; node_name[inode] = "S035B63"; node_title[inode] = "G(h- h- h+ >=0 neutrals nu(tau) (ex. K(S)0 --> pi- pi+) (``3-prong'')) / G(total)"; node_gamma[inode] = "Gamma55";
  ++inode; node_name[inode] = "S035B64"; node_title[inode] = "G(h- h- h+ nu(tau) (ex. K0)) / G(h- h- h+ >=0 neutrals nu(tau) (ex. K(S)0 --> pi- pi+) (``3-prong''))"; node_gamma[inode] = "Gamma57by55";
  ++inode; node_name[inode] = "S035B67"; node_title[inode] = "G(Kbar0 h- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma34";
  ++inode; node_name[inode] = "S035B68"; node_title[inode] = "G(Kbar0 h- pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma39";
  ++inode; node_name[inode] = "S035B69"; node_title[inode] = "G(pi- K(S)0 K(S)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma47";
  ++inode; node_name[inode] = "S035B71"; node_title[inode] = "G(h- h- h+ nu(tau) (ex. K0, omega)) / G(total)"; node_gamma[inode] = "Gamma58";
  ++inode; node_name[inode] = "S035B72"; node_title[inode] = "G(h- h- h+ 2pi0 nu(tau) (ex. K0, omega, eta)) / G(total)"; node_gamma[inode] = "Gamma77";
  ++inode; node_name[inode] = "S035B73"; node_title[inode] = "G(h- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma8";
  ++inode; node_name[inode] = "S035B74"; node_title[inode] = "G(h- 2pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma18";
  ++inode; node_name[inode] = "S035B75"; node_title[inode] = "G((particles)- >=0 neutrals >=0 K0 nu(tau) (``1-prong'')) / G(total)"; node_gamma[inode] = "Gamma1";
  ++inode; node_name[inode] = "S035B76"; node_title[inode] = "G(h- h- h+ pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma65";
  ++inode; node_name[inode] = "S035B77"; node_title[inode] = "G(h- h- h+ 2pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma75";
  ++inode; node_name[inode] = "S035B78"; node_title[inode] = "G(h- h- h+ >=1 pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma64";
  ++inode; node_name[inode] = "S035B79"; node_title[inode] = "G(h- 4pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma29";
  ++inode; node_name[inode] = "S035B89"; node_title[inode] = "G(pi- pi- pi+ eta nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma136";
  ++inode; node_name[inode] = "S035B97"; node_title[inode] = "G(h- nu(tau)) / G(e- nubar(e) nu(tau))"; node_gamma[inode] = "Gamma8by5";
  ++inode; node_name[inode] = "S035B98"; node_title[inode] = "G(Kbar0 pi- 2pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma44";
  ++inode; node_name[inode] = "S035C01"; node_title[inode] = "G(h- >= 1pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma12";
  ++inode; node_name[inode] = "S035C02"; node_title[inode] = "G(h- >= 3pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma25";
  ++inode; node_name[inode] = "S035C03"; node_title[inode] = "G(h- h- h+ >= 2pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma74";
  ++inode; node_name[inode] = "S035C1 "; node_title[inode] = "G(pi- K(S)0 K(L)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma48";
  ++inode; node_name[inode] = "S035C18"; node_title[inode] = "G(pi- pi- pi+ nu(tau)) / G(total)"; node_gamma[inode] = "Gamma59";
  ++inode; node_name[inode] = "S035C19"; node_title[inode] = "G(pi- pi- pi+ nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma60";
  ++inode; node_name[inode] = "S035C20"; node_title[inode] = "G(pi- pi- pi+ nu(tau) (ex. K0, omega)) / G(total)"; node_gamma[inode] = "Gamma62";
  ++inode; node_name[inode] = "S035C21"; node_title[inode] = "G(K- pi- pi+ nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma85";
  ++inode; node_name[inode] = "S035C22"; node_title[inode] = "G(pi- pi- pi+ pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma68";
  ++inode; node_name[inode] = "S035C23"; node_title[inode] = "G(pi- pi- pi+ pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma69";
  ++inode; node_name[inode] = "S035C24"; node_title[inode] = "G(pi- pi- pi+ pi0 nu(tau) (ex. K0, omega)) / G(total)"; node_gamma[inode] = "Gamma70";
  ++inode; node_name[inode] = "S035C25"; node_title[inode] = "G(K- pi- pi+ pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma88";
  ++inode; node_name[inode] = "S035C27"; node_title[inode] = "G(K- pi0 eta nu(tau)) / G(total)"; node_gamma[inode] = "Gamma130";
  ++inode; node_name[inode] = "S035C28"; node_title[inode] = "G(Kbar0 pi- eta nu(tau)) / G(total)"; node_gamma[inode] = "Gamma132";
  ++inode; node_name[inode] = "S035C31"; node_title[inode] = "G(K- pi- h+ nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma80";
  ++inode; node_name[inode] = "S035C32"; node_title[inode] = "G(K- pi- h+ nu(tau) (ex. K0)) / G(pi- pi- pi+ nu(tau) (ex. K0))"; node_gamma[inode] = "Gamma80by60";
  ++inode; node_name[inode] = "S035C33"; node_title[inode] = "G(K- pi- h+ pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma81";
  ++inode; node_name[inode] = "S035C34"; node_title[inode] = "G(K- pi- h+ pi0 nu(tau) (ex. K0)) / G(pi- pi- pi+ pi0 nu(tau) (ex. K0))"; node_gamma[inode] = "Gamma81by69";
  ++inode; node_name[inode] = "S035C35"; node_title[inode] = "G(pi- K- K+ nu(tau)) / G(pi- pi- pi+ nu(tau) (ex. K0))"; node_gamma[inode] = "Gamma93by60";
  ++inode; node_name[inode] = "S035C36"; node_title[inode] = "G(pi- K- K+ pi0 nu(tau)) / G(pi- pi- pi+ pi0 nu(tau) (ex. K0))"; node_gamma[inode] = "Gamma94by69";
  ++inode; node_name[inode] = "S035C37"; node_title[inode] = "G(Kbar0 pi- >=1 pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma43";
  ++inode; node_name[inode] = "S035C38"; node_title[inode] = "G(K- K0 >=0 pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma38";
  ++inode; node_name[inode] = "S035C40"; node_title[inode] = "G(K- pi- pi+ >=0 pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma83";
  ++inode; node_name[inode] = "S035C44"; node_title[inode] = "G(pi- pi0 K0 Kbar0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma49";
  ++inode; node_name[inode] = "S035C47"; node_title[inode] = "G(strange) / G(total)"; node_gamma[inode] = "Gamma110";
  ++inode; node_name[inode] = "S035C5 "; node_title[inode] = "G(Kbar0 h- h- h+ nu(tau)) / G(total)"; node_gamma[inode] = "Gamma53";
  ++inode; node_name[inode] = "S035C54"; node_title[inode] = "G(K- pi- pi+ pi0 nu(tau) (ex. K0, eta)) / G(total)"; node_gamma[inode] = "Gamma89";
  ++inode; node_name[inode] = "S035C6 "; node_title[inode] = "G(K- pi- pi+ nu(tau)) / G(total)"; node_gamma[inode] = "Gamma84";
  ++inode; node_name[inode] = "S035C61"; node_title[inode] = "G(K- omega nu(tau)) / G(total)"; node_gamma[inode] = "Gamma151";
  ++inode; node_name[inode] = "S035C7 "; node_title[inode] = "G(K- pi- pi+ pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma87";
  ++inode; node_name[inode] = "S035C8 "; node_title[inode] = "G(pi- K- K+ pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma94";
  ++inode; node_name[inode] = "S035C9 "; node_title[inode] = "G(K- K- K+ nu(tau)) / G(total)"; node_gamma[inode] = "Gamma96";
  ++inode; node_name[inode] = "S035R1 "; node_title[inode] = "G(mu- nubar(mu) nu(tau)) / G(total)"; node_gamma[inode] = "Gamma3";
  ++inode; node_name[inode] = "S035R14"; node_title[inode] = "G(h- omega nu(tau)) / G(h- h- h+ pi0 nu(tau) (ex. K0))"; node_gamma[inode] = "Gamma150by66";
  ++inode; node_name[inode] = "S035R15"; node_title[inode] = "G(h- omega >=0 neutrals nu(tau)) / G(total)"; node_gamma[inode] = "Gamma149";
  ++inode; node_name[inode] = "S035R2 "; node_title[inode] = "G(e- nubar(e) nu(tau)) / G(total)"; node_gamma[inode] = "Gamma5";
  ++inode; node_name[inode] = "S035R20"; node_title[inode] = "G(h- 2pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma19";
  ++inode; node_name[inode] = "S035R21"; node_title[inode] = "G(h- 3pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma26";
  ++inode; node_name[inode] = "S035R23"; node_title[inode] = "G(h- omega nu(tau)) / G(total)"; node_gamma[inode] = "Gamma150";
  ++inode; node_name[inode] = "S035R24"; node_title[inode] = "G((particles)- >=0 neutrals >=0 K(L)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma2";
  ++inode; node_name[inode] = "S035R26"; node_title[inode] = "G(K- >=0 pi0 >=0 K0 >=0 gamma nu(tau)) / G(total)"; node_gamma[inode] = "Gamma31";
  ++inode; node_name[inode] = "S035R27"; node_title[inode] = "G(K- >=1 (pi0 or K0 or gamma) nu(tau)) / G(total)"; node_gamma[inode] = "Gamma32";
  ++inode; node_name[inode] = "S035R28"; node_title[inode] = "G(h- h- h+ nu(tau)) / G(total)"; node_gamma[inode] = "Gamma56";
  ++inode; node_name[inode] = "S035R30"; node_title[inode] = "G(h- h- h+ >=1 neutrals nu(tau)) / G(total)"; node_gamma[inode] = "Gamma63";
  ++inode; node_name[inode] = "S035R31"; node_title[inode] = "G(h- h- h+ >=0 neutrals >=0 K(L)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma54";
  ++inode; node_name[inode] = "S035R32"; node_title[inode] = "G(pi- pi0 eta nu(tau)) / G(total)"; node_gamma[inode] = "Gamma126";
  ++inode; node_name[inode] = "S035R33"; node_title[inode] = "G(3h- 2h+ >=0 neutrals nu(tau) (ex. K(S)0 --> pi- pi+) (``5-prong'')) / G(total)"; node_gamma[inode] = "Gamma102";
  ++inode; node_name[inode] = "S035R34"; node_title[inode] = "G(K- h- h+ >=0 neutrals nu(tau)) / G(total)"; node_gamma[inode] = "Gamma79";
  ++inode; node_name[inode] = "S035R38"; node_title[inode] = "G(3h- 2h+ nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma103";
  ++inode; node_name[inode] = "S035R39"; node_title[inode] = "G(3h- 2h+ pi0 nu(tau) (ex. K0)) / G(total)"; node_gamma[inode] = "Gamma104";
  ++inode; node_name[inode] = "S035R40"; node_title[inode] = "G(pi- K- K+ nu(tau)) / G(total)"; node_gamma[inode] = "Gamma93";
  ++inode; node_name[inode] = "S035R41"; node_title[inode] = "G(K- pi- pi+ >=0 neutrals nu(tau)) / G(total)"; node_gamma[inode] = "Gamma82";
  ++inode; node_name[inode] = "S035R42"; node_title[inode] = "G(h- >=1 neutrals nu(tau)) / G(total)"; node_gamma[inode] = "Gamma11";
  ++inode; node_name[inode] = "S035R43"; node_title[inode] = "G(h- >=0 K(L)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma7";
  ++inode; node_name[inode] = "S035R44"; node_title[inode] = "G(h- >=2 pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma17";
  ++inode; node_name[inode] = "S035R46"; node_title[inode] = "G(K- K0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma37";
  ++inode; node_name[inode] = "S035R5 "; node_title[inode] = "G(mu- nubar(mu) nu(tau)) / G(e- nubar(e) nu(tau))"; node_gamma[inode] = "Gamma3by5";
  ++inode; node_name[inode] = "S035R6 "; node_title[inode] = "G(pi- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma9";
  ++inode; node_name[inode] = "S035R7 "; node_title[inode] = "G(K- nu(tau)) / G(total)"; node_gamma[inode] = "Gamma10";
  ++inode; node_name[inode] = "S035R8 "; node_title[inode] = "G(pi- pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma14";
  ++inode; node_name[inode] = "S035R84"; node_title[inode] = "G(h- pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma13";
  ++inode; node_name[inode] = "S035R97"; node_title[inode] = "G(h- >= 3pi0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma24";
  ++inode; node_name[inode] = "S035S01"; node_title[inode] = "G(total)"; node_gamma[inode] = "GammaAll";
  ++inode; node_name[inode] = "S035Y01"; node_title[inode] = "G(pi- nu(tau)) / G(e- nubar(e) nu(tau))"; node_gamma[inode] = "Gamma9by5";
  ++inode; node_name[inode] = "S035Y02"; node_title[inode] = "G(K- nu(tau)) / G(e- nubar(e) nu(tau))"; node_gamma[inode] = "Gamma10by5";
  ++inode; node_name[inode] = "S035Y03"; node_title[inode] = "G(K- nu(tau)) / G(pi- nu(tau))"; node_gamma[inode] = "Gamma10by9";
  ++inode; node_name[inode] = "S035Z01"; node_title[inode] = "G(pi- omega nu(tau)) / G(total)"; node_gamma[inode] = "Gamma800";
  ++inode; node_name[inode] = "S035Z02"; node_title[inode] = "G(K- phi nu(tau) (phi->KK)) / G(total)"; node_gamma[inode] = "Gamma801";
  ++inode; node_name[inode] = "S035Z03"; node_title[inode] = "G(K- pi- pi+ nu(tau) (ex. K0, omega)) / G(total)"; node_gamma[inode] = "Gamma802";
  ++inode; node_name[inode] = "S035Z04"; node_title[inode] = "G(K- pi- pi+ pi0 nu(tau) (ex. K0, omega, eta)) / G(total)"; node_gamma[inode] = "Gamma803";
  ++inode; node_name[inode] = "S035Z05"; node_title[inode] = "G(pi- K(L)0 K(L)0 nu(tau)) / G(total)"; node_gamma[inode] = "Gamma804";
  ++inode; node_name[inode] = "S035Z06"; node_title[inode] = "G(a1- (-> pi- gamma) nu(tau)) / G(total)"; node_gamma[inode] = "Gamma805";
  //
  int nnode=inode+1;
  for (inode=0;inode<nnode;++inode) {
    node_gamma[inode].erase(0,5);
  }
  //
  const int istart=135;
  //
  for (imeas=0;imeas<nmeas;++imeas) {
    inode=-1;
    for (inode=0;inode<nnode;++inode) {
      if (strcmp(meas_gamma[imeas].data(),node_gamma[inode].data())==0) break;
    }
    if (inode==-1 || inode==nnode) {cout << "check imeas = " << imeas << " inode = " << inode << endl; exit(1);}
    //
    cout << "* " << Form("%4d",istart+imeas) << " " << Form("%-5s",meas_gamma[imeas].data()) << "   " << Form("%7s",node_name[inode].data())
	 << Form("%12.9f",meas_val_adj[imeas]) << " " << Form("%12.6e",meas_err_adj[imeas]) << " " 
	 << Form("%6s",meas_exp[imeas].data()) << " " << Form(" AAAAAA") << "      " << "XXY" << "  " << node_title[inode] << endl;
  }
  //
  for (imeas=0;imeas<nmeas;++imeas) {
    for (jmeas=0;jmeas<nmeas;++jmeas) {
      if (jmeas>imeas) {
	if (meas_corr[imeas][jmeas]!=0) {
	  cout << "% " << Form("%4d",istart+imeas) << " " << Form("%4d",istart+jmeas) << " " << Form("%20.10f",meas_corr[imeas][jmeas]*100) << endl;
	}
      }
    }
  }
  //
  return 0;
  //
}
