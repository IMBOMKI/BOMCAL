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
  
  ///////////////////
  // Set Variables //
  ///////////////////

  /* branching ratios */
  Double_t SIGBr=1.7*TMath::Power(10,-12); // Set to current limit on mu->e+ conversion
  Double_t RMCBr=6.82*TMath::Power(10,-7);

  /* common variables */
  Double_t eff=0.01;
  Double_t protonNum=TMath::Power(10,21);
  Double_t muonStopRate=TMath::Power(10,-3);

  /* SIG specific variables */
  Double_t capturingRate=0.61; // capturing rate of Aluminum

  /* RMC specific vairables */
  Double_t PairCreationProb=0.97;
  Double_t VinTProb=0.0058; // Probability that vertex is in target
  Double_t ECut=0.016;      // Probability that positron has energy hihgher than 90.30

  /////////////////////////
  // Set Composite Model //
  /////////////////////////

  // SIG component
  // Landau distribuition determined by checking smeard energy from muon stopping target

  RooRealVar landauMean("landauMean","mean of Landau",0.000447);
  RooRealVar landauVariance("landauVariance","variance of Landau", 0.000112);
  RooLandau sig("sig", "sig p.d.f.", landauMean, landuVariance);

  Double_t sigFrac=protonNum*muonStopRate*capturingRate*SIGBr;

  // RMC component

  std::string rootPath="/home/bomki/ICEDUST/BOMKI_analysis/v999/rootfile/";
  std::string rootFile="extrmc_1e7.root";

  Double_t rmcFrac=protonNum*muonStopRate*RMCBr*PairCreationProb*VinTProb*ECut;

  TFile *f = TFile::Open(TString(rootPath+rootFile));
  TTree *t = (TTree*)f->Get("trdata");

}