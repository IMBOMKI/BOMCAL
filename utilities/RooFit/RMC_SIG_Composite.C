#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLandau.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <iostream>
#include <sstream>
#include <TMath.h>

void RMC_SIG_Composite(){
  
  // Set Variables

  /* branching ratios */
  Double_t SIGBr=1.7*TMath::Power(10,-12); // Set to current limit on mu->e+ conversion
  Double_t RMCBr=6.82*TMath::Power(10,-7);

  /* common variables */
  Double_t eff=0.01;
  Double_t protonNum=TMath::Power(10,21);
  Double_t muStopRate=TMath::Power(10,-3);

  /* RMC specific vairables */
  Double_t PairCreationProb=0.97;


  // SIGNAL component

  RooRealVar landauMean("landauMean","mean of Landau",0.000447);
  RooRealVar landauVariance("landauVariance","variance of Landau", 0.000112);
  
  RooRealVar sigFrac("sigFrac", "fraction of signal", 


  // RMC component

  std::string rootPath="/home/bomki/ICEDUST/BOMKI_analysis/v999/rootfile/";
  std::string rootFile="extrmc_1e7.root";

  TFile *f = TFile::Open(TString(rootPath+rootFile));
  TTree *t = (TTree*)f->Get("trdata");

}
