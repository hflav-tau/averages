//
// aluAverPlot.cc
// root -b -q "aluAverPlot.cc+(\"plot.input\")"
//
// plot measurements and averages listed in ASCII .input file
// inspired by S.Banerjee code for HFAG-tau plots
//
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
#include "TPaveText.h"

#include "HFAGTauLabel.cc"

void PlotStyleSetup(Float_t scale = 1, Int_t canvas_width = 560)
{
  gROOT->SetStyle("Pub");

  //--- canvas size
  gStyle->SetCanvasDefW(canvas_width);
  //--- height will be recalculated later
  gStyle->SetCanvasDefH(Int_t(0.75*gStyle->GetCanvasDefW()));
  
  //--- normal plots
  // gStyle->SetPadLeftMargin(0.20);
  // gStyle->SetPadRightMargin(0.06);
  // gStyle->SetPadTopMargin(0.08);
  // gStyle->SetPadBottomMargin(0.15);

  //--- average plots
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.50);

  gStyle->SetLineWidth(3);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetEndErrorSize(6);

  gStyle->SetFuncColor(kBlue);
  gStyle->SetLineColor(kBlack);
  gStyle->SetStatColor(kWhite);

  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetTitleOffset(1.6, "y");

  gStyle->SetLabelOffset(0.01, "x");
  gStyle->SetTickLength(0.06*scale, "x");
  gStyle->SetNdivisions(505, "x");

  //--- do not print y axis
  gStyle->SetAxisColor(kWhite, "y");
  gStyle->SetLabelSize(0, "y");
  gStyle->SetLabelColor(kWhite, "y");
  gStyle->SetTickLength(0, "y");

  gStyle->SetTextSize(0.04*scale);
  gStyle->SetTitleSize(0.05*scale, "x");
  gStyle->SetLabelSize(0.04*scale, "x");

  const Int_t plotFont(42);
  gStyle->SetLabelFont(plotFont, "xyz");
  gStyle->SetTitleFont(plotFont, "xyz");
  gStyle->SetTextFont(plotFont);
  gStyle->SetStatFont(plotFont);

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.5);
  gStyle->SetMarkerColor(kBlack);

  gStyle->SetStatH(0.10);
  gStyle->SetStatW(0.16);
  gStyle->SetFitFormat(".4g");
  gStyle->SetStatFormat(".4g");

  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(000000);

  //--- average specific
  gStyle->SetTitleTextColor(kBlack);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetCanvasBorderSize(0);

  gROOT->ForceStyle();
}

// ////////////////////////////////////////////////////////////////////////////

void aluAverPlot(const std::string& filename = "plot.input",
		 const Int_t nPoints_def = 6,
		 Int_t canvas_width = 560);

void aluAverPlot(const std::string& filename, const Int_t nPoints_def, Int_t canvas_width)
{  
  //--- num of points in plot
  Int_t nPoints(0);
  const Int_t nPointsMax(20);

  //--- fixed vertical offset
  const Float_t fnPoints_offset(0.6);
  
  TString title, expname[nPointsMax], tempstring;
  TString fitLabel("Summer 2014");
  Float_t precision;

  // -- set up frame for plot
  Float_t fXmin;
  Float_t fXmax;
  Float_t fYmin(0);
  Float_t fYmax;
  
  Float_t meas[nPointsMax];
  Float_t stath[nPointsMax], statl[nPointsMax], stat[nPointsMax];
  Float_t systh[nPointsMax], systl[nPointsMax], syst[nPointsMax];
  Float_t errorh[nPointsMax], errorl[nPointsMax], error[nPointsMax];
  Float_t meascl[nPointsMax];
  Float_t sfact[nPointsMax];

  TString measinfo[nPointsMax];
  
  //--- open input file
  ifstream ifs(filename.c_str()) ; if (!ifs.good()) {cout << "Cannot open input file '" << filename << "'" << endl ; exit(1) ;}
  char buffer[200] ;
  Int_t ref_meas(0);
  
  while(ifs.good()) {
    if (ifs.eof()) break;
    char first_ch; ifs.get(first_ch) ;

    if (first_ch=='#' || first_ch=='\n') { // Skip this line
      ifs.ignore(200, '\n') ;

    } else if (first_ch=='i') {
      //--- plot limits and title: i <xmin> <xmax> <format> <title>
      ifs >> fXmin >> fXmax >> precision;
      ifs.getline(buffer, 200, '\n');
      title = TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace
    } else if (first_ch=='f') {
      //--- fit label
      ifs.getline(buffer, 200, '\n');
      fitLabel = TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace
    } else if (first_ch=='h' || first_ch=='p') {
      //--- other average: h <mean> <sigma> <s-factor> <CL> <label>
      //--- pdg average: p <mean> <sigma> <s-factor> <CL> <label> -- other average
      if (first_ch=='p') ref_meas = nPoints;
      //--- get averages: mean, sigma, sfact, CL
      measinfo[nPoints] = first_ch;
      ifs >>  meas[nPoints] >> stath[nPoints] >> sfact[nPoints] >> meascl[nPoints];
      statl[nPoints] = stath[nPoints];
      systh[nPoints] = systl[nPoints] = 0;
      ifs.getline(buffer,200,'\n');
      expname[nPoints++] = TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace

    } else if (first_ch=='m') {
      //
      // measurements
      // - m   symm errors
      // - ma  asymm errors
      // - mt  "this work" symm errors
      // - mat "this work" asymm errors
      //
      measinfo[nPoints] = first_ch;
      bool asymm_input = false;
      char ch2 ; ifs.get(ch2);

      if (ch2 == 'a') {
	asymm_input = true;
	char ch3 ; ifs.get(ch3);
	if (ch3=='t') {
	  measinfo[nPoints] = 't';
	} else {
	  ifs.putback(ch3);
	}
      } else if (ch2 == 't') {
	measinfo[nPoints] = 't';
      } else {
	ifs.putback(ch2);
      }

      //--- parse input
      if (asymm_input) {
	ifs >> meas[nPoints] >> stath[nPoints] >> statl[nPoints] >> systh[nPoints] >> systl[nPoints];
	stath[nPoints] = fabs(stath[nPoints]);
	statl[nPoints] = fabs(statl[nPoints]);
	systh[nPoints] = fabs(systh[nPoints]);
	systl[nPoints] = fabs(systl[nPoints]);
      } else {
	ifs >> meas[nPoints] >> stath[nPoints] >> systh[nPoints];
	statl[nPoints] = stath[nPoints];
	systl[nPoints] = systh[nPoints];
      }

      ifs.getline(buffer,200,'\n') ;
      expname[nPoints] = TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace
      if (expname[nPoints].Length()) nPoints++;

    } else { // discard
      ifs.ignore(200,'\n') ;
    }
  }
  cout << "read " << nPoints << " data lines from file '" << filename << "'" << endl ;
  
  Int_t nexp(0);
  Float_t num(0), den(0), aver(0), err(0), wt;
  Float_t meas_xmin(meas[0]);
  Float_t meas_xmax(meas[0]);
  Float_t meas_width_min(stath[0]+statl[0]+systh[0]+systl[0]);
  Float_t meas_width_xmin(meas[0]);
  Float_t meas_width_xmax(meas[0]);

  for (Int_t i=0;i<nPoints;++i) {

    stat[i] = 0.5*(stath[i]+statl[i]);
    syst[i] = 0.5*(systh[i]+systl[i]);

    errorh[i] = sqrt(pow(stath[i], 2) + pow(systh[i], 2));
    errorl[i] = sqrt(pow(statl[i], 2) + pow(systl[i], 2));
    error[i] = 0.5*(errorh[i]+errorl[i]);

    if (measinfo[i] == "h" || measinfo[i] == "p") {
      if (measinfo[i] == "h") {
	cout << i << " hfag ave  ";
      } else {
	cout << i << " pdg  ave  ";
      }
      cout << meas[i] << " +- " << stath[i] << " sfact " << sfact[i] << " CL " << meascl[i]
	   << " " << expname[i] << endl;
    } else if (measinfo[i] == "m" || measinfo[i] == "t") {
      if (measinfo[i] == "m") {
	cout << i << " meas      " << meas[i];
      } else {
	cout << i << " this work " << meas[i];
      }
      if (stath[i] == statl[i]) {
	cout << " +- " << stath[i];
      } else {
	cout << " +" << stath[i] << "-" << statl[i];
      }
      if (systh[i] == systl[i]) {
	cout << " +- " << systh[i];
      } else {
	cout << " +" << systh[i] << "-" << systl[i];
      }
      cout << " " << expname[i] << endl;

      wt=1./(stat[i]*stat[i]+syst[i]*syst[i]);
      num+=meas[i]*wt;
      den+=wt;
      nexp++;
    } else {
      cout << "+++ error, expname = " << expname[i] << endl;
    }

    //--- compute min/max x position of drawn bars
    Float_t xmintmp = meas[i] - errorl[i];
    Float_t xmaxtmp = meas[i] + errorh[i];
    if (xmintmp < meas_xmin) meas_xmin = xmintmp;
    if (xmaxtmp > meas_xmax) meas_xmax = xmaxtmp;

    Float_t meas_width_min_tmp(errorh[i]+errorl[i]);
    // cout <<  meas_width_min_tmp << " " << meas_width_min << endl;
    if (meas_width_min_tmp < meas_width_min) {
      meas_width_min = meas_width_min_tmp;
      meas_width_xmin = meas[i] - errorl[i];
      meas_width_xmax = meas[i] + errorh[i];
    }
  }

  if (den>0) {
    aver=num/den;
    err =sqrt(1/den);
  } else {
    aver = 0;
    err = 1;
  }

  //--- if xmin >= xmax then use auto-computed x axis limits
  if (fXmin >= fXmax) {
    //--- set x axis limits where drawing arrives plus 10% of the total width
    fXmin = meas_xmin - 0.10 * (meas_xmax - meas_xmin);
    fXmax = meas_xmax + 0.10 * (meas_xmax - meas_xmin);

    //--- avoid excessive x axis limits when errors are large
    const Float_t max_width_mult(8);
    if (fXmin < meas_width_xmin - max_width_mult*meas_width_min) {
      fXmin = meas_width_xmin - max_width_mult*meas_width_min;
    }
    if (fXmax > meas_width_xmax + max_width_mult*meas_width_min) {
      fXmax = meas_width_xmax + max_width_mult*meas_width_min;
    }
    // cout << meas_width_min << " " << meas_width_xmax << " " << meas_width_xmax + max_width_mult*meas_width_min << endl;
  }

  //--- compute plot ymax according to number to input points
  fYmax = Float_t(max(nPoints, nPoints_def)) + fnPoints_offset;
  TString sprecision = Form("%s%3.1f%s", "%", precision, "f");

  cout << "x axis min/max:       " << fXmin << " " << fXmax << endl
       << "measurements min/max: " << meas_xmax << " " << meas_xmin << endl
       << "total num of points:  " << nPoints << " ymin/ymax: " << fYmin << " " << fYmax << endl
       << "measurements num:     " << nexp << ", average: " << aver << " +- " << err << endl
       << "format precision:     " << sprecision << endl
    ;

  //
  // Canvas
  //
  PlotStyleSetup(1, canvas_width);

  Float_t plotScaleFact;
  TCanvas* canvas;
  {
    //--- vertical slot in pixel needed for each shown point
    const Float_t vert_slot(0.08 * gStyle->GetCanvasDefW());

    //--- float default number of points
    const Float_t fnPoints_def(nPoints_def);
    const Float_t fnPoints(nPoints);
    
    const Float_t padTopMarginAbs(0.3 * vert_slot);
    const Float_t padBottomMarginAbs( (title == "" ? 0.5625 : 1.5 ) * vert_slot );

    const Float_t plot_height(vert_slot * (max(fnPoints, fnPoints_def) + fnPoints_offset) );

    const Int_t canvas_height( TMath::Nint(plot_height + padTopMarginAbs + padBottomMarginAbs) );

    plotScaleFact =
      0.75 * Float_t(gStyle->GetCanvasDefW()) /
      Float_t(min( gStyle->GetCanvasDefW(), canvas_height));

    // cout << "plotScaleFact " << plotScaleFact << endl;
    // cout << " " <<  gStyle->GetCanvasDefH() << " " << canvas_height << " " << vert_slot << " " << plot_height << endl;

    PlotStyleSetup(plotScaleFact, canvas_width);

    gStyle->SetCanvasDefH(canvas_height);
    gStyle->SetPadTopMargin(padTopMarginAbs / Float_t(canvas_height));
    gStyle->SetPadBottomMargin(padBottomMarginAbs / Float_t(canvas_height));

    canvas = new TCanvas("canvas", "canvas");
  }
  // canvas->SetGrid(1, 1);

  cout << "Drawing frame: " << fXmin << ", " << fYmin << ", " << fXmax << ", " << fYmax
       << ", x axis title=\"" << title << "\"" << endl;

  TH1F *hpx = canvas->DrawFrame(fXmin, fYmin, fXmax, fYmax, "");
  hpx->GetXaxis()->CenterTitle(kTRUE);
  hpx->SetXTitle(title);

  //--- draw band for average measurement
  TBox average_band;
  average_band.SetFillStyle(1000);
  // average_band.SetFillColor(kGreen-9);
  // average_band.SetFillColor(kOrange-9);
  average_band.SetFillColor(kYellow-9);
  average_band.DrawBox(meas[ref_meas]-statl[ref_meas], fYmin,
		       meas[ref_meas]+stath[ref_meas], fYmax);

  //--- redraw axis above the average band
  canvas->RedrawAxis();

  Float_t y[nPoints], ey[nPoints], eyl[nPoints], eyh[nPoints];
  Int_t fColor[nPoints], fSymbol[nPoints];

  for (Int_t i=0; i<nPoints; ++i) {
    //--- vertical position
    y[i] = Float_t(nPoints-i) -1 + 0.8;
    ey[i] = 0.0;
    eyl[i] = 0.0;
    eyh[i] = 0.0;
    fColor[i]  = kBlack;
    fSymbol[i] = 20;
  }

  //--- x label position at 107% of plot width
  const Float_t fTxtX(0.07);
  Float_t xtext(fXmax + fTxtX*(fXmax-fXmin));

  //--- y label base position at position dependent on plot height
  Float_t ytext;

  //--- result label
  TLatex  tlexp;
  tlexp.SetTextSize(gStyle->GetLabelSize());

  //--- result value
  TLatex  tlval;
  tlval.SetTextSize(gStyle->GetLabelSize()*1.00);

  Int_t color;
  for (Int_t i=0; i<nPoints; ++i) {
    color = kBlack;
    if (measinfo[i] == "p") {
      color = kBlue;
    } else if (measinfo[i] == "t") {
      color = kRed;
    }

    //--- total errors
    TGraph* GraphA = new TGraphAsymmErrors(1, &meas[i], &y[i], &errorl[i], &errorh[i], &eyl[i], &eyh[i]);
    GraphA->SetLineColor(color);
    GraphA->SetMarkerColor(color);
    GraphA->SetMarkerStyle(fSymbol[i]);
    GraphA->Draw("p");

    //--- stat errors
    if (systl[i] != 0 || systh[i] != 0) {
      TGraph* GraphB = new TGraphAsymmErrors(1, &meas[i], &y[i], &statl[i], &stath[i], &eyl[i], &eyh[i]);
      GraphB->SetLineColor(color);
      GraphB->Draw("]");
    }

    tlexp.SetTextColor(color);
    tlval.SetTextColor(color);

    ytext = y[i];
    tlexp.DrawLatex(xtext, ytext, expname[i]);
    // cout << "ytext " << ytext << " " << i << endl;

    TString valtxt(Form(sprecision.Data(),meas[i]));

    if (stath[i] == statl[i]) {
      valtxt += " #pm ";
      valtxt += Form(sprecision.Data(),stath[i]);
    } else {
      valtxt += " #splitline{+";
      valtxt += Form(sprecision.Data(),stath[i]);
      valtxt += "}{- ";
      valtxt += Form(sprecision.Data(),statl[i]);
      valtxt += "}";
    }

    if (systh[i] != 0 || systl[i] != 0) {
      if (systh[i] == systl[i]) {
	valtxt += " #pm ";
	valtxt += Form(sprecision.Data(),systh[i]);
      } else {
	valtxt += " #splitline{+";
	valtxt += Form(sprecision.Data(),systh[i]);
	valtxt += "}{- ";
	valtxt += Form(sprecision.Data(),systl[i]);
	valtxt += "}";
      }
    }

    if (measinfo[i] == "h" || measinfo[i] == "p") {
      if (sfact[i] >0) { // S-factor
	valtxt += Form(" (S = %3.1f)", sfact[i]);
      } else if (meascl[i] >0) {
	valtxt += Form(" (CL = %4.1f%%)", 100*meascl[i]);
      }
    }

    tlval.DrawLatex(xtext, ytext-0.44, valtxt);
    // cout << "value " << valtxt << endl;
  }

  const Float_t xtextNDC( (xtext-TVirtualPad::Pad()->GetX1()) / (TVirtualPad::Pad()->GetX2() - TVirtualPad::Pad()->GetX1()) );
  if (fitLabel != "") {
    HFAGTauLabel(fitLabel, xtextNDC, 0.03*plotScaleFact, 0.9 * Float_t(canvas_width) / Float_t(560));
  }
  std::string basefname(filename);
  size_t extPos = basefname.rfind('.');
  if (extPos != std::string::npos) {
    //--- Erase the current extension.
    basefname.assign(filename, 0, extPos);
  }

  canvas->SaveAs(std::string(basefname + ".pdf").c_str());
  canvas->SaveAs(std::string(basefname + ".eps").c_str());
  // canvas->SaveAs(std::string(basefname + ".gif").c_str());
  canvas->SaveAs(std::string(basefname + ".png").c_str());
  // canvas->Clear();

  return;
}
