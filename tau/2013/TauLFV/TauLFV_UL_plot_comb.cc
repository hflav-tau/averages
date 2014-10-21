///////////////////////////////////////////////////////////////////////////////
//
// TauLFV_UL_plot.cc
//
// produce tau LFV upper limits plot
// execute as follows:
//
// > root -l -b -q TauLFV_UL_plot.cc+
//
// This will create a .eps file with the plot.
//
//                          Last update  2013/07/25
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TLine.h"
#include "TBox.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "../plots/HFAGTauLabel.cc"

enum DECAY {
  EGAMMA = 0,
  MUGAMMA = 1,

  /*
  EPI0 = 2,
  MUPI0 = 3,
  EETA = 4,
  MUETA = 5,
  EETAP = 6,
  MUETAP = 7,
  */
  EKS0 = 2,
  MUKS0 = 3,
  /* 
 EF0 = 10, 
  MUF0 = 11,
  */
  ERHO = 4,   
 MURHO = 5,
  EKSTAR = 6,
  MUKSTAR = 7,
  EAKSTAR = 8,
  MUAKSTAR = 9,
  EPHI = 10,
  MUPHI = 11,
  EOMEGA = 12,
  MUOMEGA = 13,
  EEE = 14, //order of charge -+-     
  MEE = 15, 
  EMM = 16,
  MMM = 17, 
  EME = 18,   
  MEM = 19, 
  /* 

 EPIPI = 28,
  MUPIPI = 29,
  EPIK = 30,
  MUPIK = 31,
  EKPI = 32,
  MUKPI = 33,
  EKK = 34,
  MUKK = 35,
  EKS0KS0 = 36,
  MUKS0KS0 = 37,
  PIEPI = 38,
  PIMUPI = 39,
  PIEK = 40,
  PIMUK = 41,
  KEK = 42,
  KMUK = 43,
  */
 PILAM = 20,
  PILAMBAR = 21,
  KLAM = 22,
  KLAMBAR = 23,
  /*
  PMUMUOS = 48,
  PMUMUSS = 49,
  */ 
 NDECAY = 24
  
};

void fillBelle(double * array)
{
  // at 10^-8
  array[EGAMMA]  = 22.;
  array[MUGAMMA] = 17.; // added .2 for clarity in display
    
  array[EKS0]    = 1.8;//2.6;
  array[MUKS0]   = 1.7;//2.3;
  
  
  array[ERHO]    = 1.9;
  array[MURHO]   = 2.1;
  array[EKSTAR]  = 3.4;
  array[MUKSTAR] = 6.6;
  array[EAKSTAR] = 3.3;
  array[MUAKSTAR]= 6.3;
  array[EPHI]    = 3.5;
  array[MUPHI]   = 7.6;
  array[EOMEGA]  = 5.2;
  array[MUOMEGA] = 6.1;
  
  array[EEE]     = 2.7;
  array[MEE]     = 1.7;
  array[EMM]     = 2.6; 
  array[MMM]     = 2.1; 
  array[EME]     = 1.4; 
  array[MEM]     = 1.6; 
  
  array[PILAM]   = 3.2;  
  array[PILAMBAR]= 2.9;  
  array[KLAM]    = 4.4;   
  array[KLAMBAR] = 3.4;  

}


void fillBaBar(double * array)
{
  // at 10^-8
  array[EGAMMA]  = 6.1;
  array[MUGAMMA] = 5.9;
    
  array[EKS0]    = 4.7;
  array[MUKS0]   = 6.9;
    
  
  array[ERHO]    = 5.2;
  array[MURHO]   = 6.2;
  array[EKSTAR]  = 6.1;
  array[MUKSTAR] = 17.0;
  array[EAKSTAR] = 5.6;
  array[MUAKSTAR]= 9.1;
  array[EPHI]    = 4.3;
  array[MUPHI]   = 18.0;
  array[EOMEGA]  = 9.4;
  array[MUOMEGA] = 10.;
  
  array[EEE]     = 3.1;
  array[MEE]     = 3.0;
  array[EMM]     = 4.1;
  array[MMM]     = 4.0;
  array[EME]     = 2.1;
  array[MEM]     = 2.6;
  
  array[PILAM]   = 6.7;
  array[PILAMBAR]= 6.5;
  array[KLAM]    = 9.2;
  array[KLAMBAR] = 5.0;
}

void fillLHCb(double * array)
{   
  // at 10^-8
  array[MMM]     = 4.6;
  //  array[PMUMUOS]   = 33.;
  //  array[PMUMUSS]   = 44.;
}
void fillHFAG_CLs(double * array)        
{                                     
  array[EGAMMA] = 5.4;
  array[MUGAMMA] = 5.0;

  array[EKS0]    = 1.4;  
  array[MUKS0]   = 1.5;  

  
  array[ERHO] = 1.5;
  array[MURHO] = 1.5;
  array[EKSTAR]= 2.3;
  array[MUKSTAR]=6.0;
  array[EAKSTAR]=2.2;
  array[MUAKSTAR]= 4.2;
  array[EPHI]=2.0;
  array[MUPHI]=6.8;
  array[EOMEGA]  = 3.3; 
  array[MUOMEGA] = 4.0;    


  array[EEE] = 1.4;
  array[MEE]= 1.1;
  array[EMM]= 1.6;
  array[MMM] = 1.2;
  array[EME] = 0.84;
  array[MEM] = 0.98;


  array[PILAM]= 1.9;
  array[PILAMBAR]= 1.8;
  array[KLAM] = 3.7;
  array[KLAMBAR] =2.0;
   
   
}

void setLabels(TH1* hist)
{
  const Int_t label_offs(2);

  TAxis * axis = hist->GetXaxis();

  axis->SetBinLabel(label_offs - 1 ," ");
  axis->SetBinLabel(label_offs + EGAMMA ,"e^{-} #gamma");
  axis->SetBinLabel(label_offs + MUGAMMA ,"#mu^{-} #gamma");
  //  axis->SetBinLabel(label_offs + EPI0    ,"e^{-} #pi^{0}");
  // axis->SetBinLabel(label_offs + MUPI0   ,"#mu^{-} #pi^{0}");
  // axis->SetBinLabel(label_offs + EETA    ,"e^{-} #eta");
  // axis->SetBinLabel(label_offs + MUETA   ,"#mu^{-} #eta");
  //axis->SetBinLabel(label_offs + EETAP   ,"e^{-} #eta'");
  // axis->SetBinLabel(label_offs + MUETAP  ,"#mu^{-} #eta'");
  axis->SetBinLabel(label_offs + EKS0    ,"e^{-} K_{S}^{0}");
  axis->SetBinLabel(label_offs + MUKS0   ,"#mu^{-} K_{S}^{0}");
  // axis->SetBinLabel(label_offs + EF0     ,"e^{-} f_{0}");
  // axis->SetBinLabel(label_offs + MUF0    ,"#mu^{-} f_{0}");
  axis->SetBinLabel(label_offs + ERHO    ,"e^{-} #rho_{0}");
  axis->SetBinLabel(label_offs + MURHO   ,"#mu^{-} #rho_{0}");
  axis->SetBinLabel(label_offs + EKSTAR  ,"e^{-} K*");
  axis->SetBinLabel(label_offs + MUKSTAR ,"#mu^{-} K*");
  axis->SetBinLabel(label_offs + EAKSTAR ,"e^{-} #bar{K*}");
  axis->SetBinLabel(label_offs + MUAKSTAR,"#mu^{-} #bar{K*}");
  axis->SetBinLabel(label_offs + EPHI    ,"e^{-} #phi");
  axis->SetBinLabel(label_offs + MUPHI   ,"#mu^{-} #phi");
  axis->SetBinLabel(label_offs + EOMEGA  ,"e^{-} #omega");
  axis->SetBinLabel(label_offs + MUOMEGA ,"#mu^{-} #omega");
  axis->SetBinLabel(label_offs + EEE,     "e^{-} e^{+} e^{-}");
  axis->SetBinLabel(label_offs + MEE,     "#mu^{-} e^{+} e^{-}");
  axis->SetBinLabel(label_offs + EMM,     "e^{-} #mu^{+} #mu^{-}");
  axis->SetBinLabel(label_offs + MMM,     "#mu^{-} #mu^{+} #mu^{-}");
  axis->SetBinLabel(label_offs + EME,     "e^{-} #mu^{+} e^{-}");
  axis->SetBinLabel(label_offs + MEM,     "#mu^{-} e^{+} #mu^{-}");
  /*
  axis->SetBinLabel(label_offs + EPIPI   ,"e^{-} #pi^{+} #pi^{-}");
  axis->SetBinLabel(label_offs + MUPIPI  ,"#mu^{-} #pi^{+} #pi^{-}");
  axis->SetBinLabel(label_offs + EPIK    ,"e^{-} #pi^{+} K^{-}");
  axis->SetBinLabel(label_offs + MUPIK   ,"#mu^{-} #pi^{+} K^{-}");
  axis->SetBinLabel(label_offs + EKPI    ,"e^{-} K^{+} #pi^{-}");
  axis->SetBinLabel(label_offs + MUKPI   ,"#mu^{-} K^{+} #pi^{-}");
  axis->SetBinLabel(label_offs + EKK     ,"e^{-} K^{+} K^{-}");
  axis->SetBinLabel(label_offs + MUKK    ,"#mu^{-} K^{+} K^{-}");
  axis->SetBinLabel(label_offs + EKS0KS0 ,"e^{-} K_{S}^{0} K_{S}^{0}");
  axis->SetBinLabel(label_offs + MUKS0KS0,"#mu^{-} K_{S}^{0} K_{S}^{0}");
  axis->SetBinLabel(label_offs + PIEPI   ,"#pi^{-} e^{+} #pi^{-}");
  axis->SetBinLabel(label_offs + PIMUPI  ,"#pi^{-} #mu^{+} #pi^{-}");
  axis->SetBinLabel(label_offs + PIEK    ,"#pi^{-} e^{+} K^{-}");
  axis->SetBinLabel(label_offs + PIMUK   ,"#pi^{-} #mu^{+} K^{-}");
  axis->SetBinLabel(label_offs + KEK     ,"K^{-} e^{+} K^{-}");
  axis->SetBinLabel(label_offs + KMUK    ,"K^{-} #mu^{+} K^{-}");
  */
  axis->SetBinLabel(label_offs + PILAM   ,"#pi^{-} #Lambda");
  axis->SetBinLabel(label_offs + PILAMBAR,"#pi^{-} #bar{#Lambda}");
  axis->SetBinLabel(label_offs + KLAM   , "K^{-} #Lambda");
  axis->SetBinLabel(label_offs + KLAMBAR, "K^{-} #bar{#Lambda}");
  //  axis->SetBinLabel(label_offs + PMUMUOS, "#bar{p} #mu^{-} #mu^{+}");
  // axis->SetBinLabel(label_offs + PMUMUSS, "p #mu^{-} #mu^{-}"); 
}

void SetUp()
{
  // gROOT->SetBatch(kTRUE);

  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  gStyle->SetStatStyle(0);
  gStyle->SetOptStat(00000000);
  
  const int plotFont(44);
  gStyle->SetLabelFont(plotFont, "xyz");
  gStyle->SetTitleFont(plotFont, "xyz");
  gStyle->SetTextFont(42);
  gStyle->SetStatFont(42);  

  gStyle->SetTitleSize(42,"xyz");
  gStyle->SetLabelSize(42,"xyz");

  gStyle->SetTitleOffset(0.5, "x");
  gStyle->SetTitleOffset(0.8, "y");
  
  gStyle->SetLabelOffset(0.006, "x");
  gStyle->SetLabelOffset(0.001, "y");
  
  gStyle->SetTickLength(0.03, "x");
  gStyle->SetTickLength(0.02, "y");
  gStyle->SetNdivisions(5, "y");
  
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(0); //not for x-axis
  gStyle->SetPadTickY(1);
  gStyle->SetErrorX(0.001);

  gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(1);
}

void TauLFV_UL_plot_comb(Int_t when=2013001)
{
  SetUp();
  unsigned int ibin;
  
  double Belle[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) Belle[ibin] = 0;
  
  double BaBar[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) BaBar[ibin] = 0;
  
  double CLEO[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) CLEO[ibin] = 0;
  
  double LHCb[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) LHCb[ibin] = 0;
  double HFAG_CLs[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) HFAG_CLs[ibin] = 0; 

  double histbins[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) histbins[ibin] = ibin;
  

  fillBelle(Belle);
  fillBaBar(BaBar);
  fillLHCb(LHCb);
  fillHFAG_CLs(HFAG_CLs);
  
  //  fillCLEO(CLEO);
  
  TH1F * hBelle = new TH1F("hBelle", "", NDECAY+2, -1.5, float(NDECAY)+0.5);
  TH1F * hBaBar = new TH1F("hBaBar", "", NDECAY+2, -1.5, float(NDECAY)+0.5);
  //  TH1F * hCLEO  = new TH1F("hCLEO",  "", NDECAY+2, -1.5, float(NDECAY)+0.5);
  TH1F * hLHCb  = new TH1F("hLHCb",  "", NDECAY+2, -1.5, float(NDECAY)+0.5); 
  TH1F * hHFAG_CLs = new TH1F("hFAG_CLs",  "", NDECAY+2, -1.5, float(NDECAY)+0.5);

  hBelle->FillN(NDECAY, histbins, Belle);
  hBaBar->FillN(NDECAY, histbins, BaBar);
  //hCLEO-> FillN(NDECAY, histbins, CLEO);
  hLHCb-> FillN(NDECAY, histbins, LHCb);
  hHFAG_CLs-> FillN(NDECAY, histbins, HFAG_CLs);  

  hBelle->Scale(1.e-8);
  hBaBar->Scale(1.e-8);
  // hCLEO->Scale(1.e-6);
  hLHCb->Scale(1.e-8);
  hHFAG_CLs->Scale(1.e-8); 

  setLabels(hBelle);
  
  hBelle->GetYaxis()->SetTitle("90% C.L. upper limits for LFV #tau decays");
  hBelle->LabelsOption("v");
  hBelle->GetYaxis()->SetLabelSize(36);

  // hBelle->SetMarkerStyle(kOpenCircle); hBelle->SetMarkerColor(kRed); hBelle->SetMarkerSize(1.0);
  // hBaBar->SetMarkerStyle(kOpenDiamond); hBaBar->SetMarkerColor(kBlue); hBaBar->SetMarkerSize(1.3);
  // hCLEO->SetMarkerStyle(kOpenSquare);  hCLEO->SetMarkerColor(kGreen);  hCLEO->SetMarkerSize(0.8);
  
  hBelle->SetMarkerStyle(kFullTriangleUp); hBelle->SetMarkerColor(kMagenta); hBelle->SetMarkerSize(1.5);
  hBaBar->SetMarkerStyle(kFullTriangleDown); hBaBar->SetMarkerColor(kCyan); hBaBar->SetMarkerSize(1.5);
  //hCLEO->SetMarkerStyle(kFullCircle);  hCLEO->SetMarkerColor(kMagenta);  hCLEO->SetMarkerSize(1.5);
  hLHCb->SetMarkerStyle(kFullSquare); hLHCb->SetMarkerColor(kBlack); hLHCb->SetMarkerSize(1.2); 
  hHFAG_CLs->SetMarkerStyle(34); hHFAG_CLs->SetMarkerColor(kGreen); hHFAG_CLs-> SetMarkerSize(1.5); 

  TCanvas *c1 = new TCanvas("c1","",1700,800); 
  c1->SetBottomMargin(0.17);
  c1->SetTopMargin(0.02);
  c1->SetLeftMargin(0.10);
  c1->SetRightMargin(0.15);
  c1->SetLogy(1) ;
  
  double ul_min=4.e-5;
  hBelle->SetMaximum(ul_min);
  double ul_max=4.e-9;
  hBelle->SetMinimum(ul_max);
  hBelle->Draw("p");
  c1->SetGridy();
  hBaBar->Draw("p,same");
  // hCLEO->Draw("p,same");
  hLHCb->Draw("p,same"); 
  hHFAG_CLs->Draw("p,same");
  c1->Update();
  
  double y_latex=2.e-5;
  
  TLine l0(EGAMMA-0.5, ul_min, EGAMMA-0.5, ul_max); l0.SetLineColor(kGray); l0.Draw();
  

  
  TLine l1(MUGAMMA+0.5,ul_min,MUGAMMA+0.5,ul_max); l1.SetLineColor(kGray); l1.Draw();
  TLatex t1(((MUGAMMA-EGAMMA)*1.0/2.0)+EGAMMA*1.0,y_latex,"l#gamma"); t1.SetTextAlign(21); t1.SetTextFont(42); t1.Draw();
  


  
  TLine l4(MUOMEGA+0.5,ul_min,MUOMEGA+0.5,ul_max); l4.SetLineColor(kGray); l4.Draw();
  TLatex t4(((MUOMEGA-EKS0 )*1.0/2.0)+EKS0*1.0,y_latex,"lh"); t4.SetTextAlign(21); t4.SetTextFont(42); t4.Draw();
  
  TLine l5(MEM+0.5,ul_min,MEM+0.5,ul_max); l5.SetLineColor(kGray); l5.Draw();
  TLatex t5(((MEM-EEE)*1.0/2.0)+EEE*1.0,y_latex,"lll"); t5.SetTextAlign(21); t5.SetTextFont(42); t5.Draw();
  

  TLine l7(KLAMBAR+0.5,ul_min,KLAMBAR+0.5,ul_max); l7.SetLineColor(kGray); l7.Draw();
  TLatex t7(((PILAM-KLAMBAR)*1.0/2.0)+KLAMBAR*1.0,y_latex,"h #Lambda"); t7.SetTextAlign(21); t7.SetTextFont(42); t7.Draw();

  TLegend *leg = new TLegend(0.85,0.3,0.95,0.5);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->SetMargin(0.5);
  //  leg->AddEntry(hCLEO,"CLEO","p");
  leg->AddEntry(hBaBar,"BaBar","p");
  leg->AddEntry(hBelle,"Belle","p");
  leg->AddEntry(hLHCb,"LHCb","p");
  leg->AddEntry(hHFAG_CLs, "HFAG", "p");
  leg->Draw();
  c1->Update();
  //
  HFAGTauLabel("2014, prelim.", -0.01, -0.02, 1.3);
  //    HFAGTauLabel("Summer 2013", -0.01, -0.02, 1.3);
    c1->Update();  
  //  HFAGTauLabel("Winter 2012", -0.01, -0.02, 1.3);
  //  c1->Update();  
  
  c1->SaveAs(Form("TauLFV_UL_%d_averaged.pdf",when));
  c1->SaveAs(Form("TauLFV_UL_%d_averaged.eps",when));
  c1->SaveAs(Form("TauLFV_UL_%d_averaged.png",when));
}
