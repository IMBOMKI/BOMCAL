#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLandau.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooBinning.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "RooTFnBinding.h"
#include "RooPlot.h"
#include "Math/DistFunc.h"
#include <iostream>
#include <sstream>
#include <TMath.h>

using namespace RooFit;

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
  Double_t lowerBound = 90.30;
  Double_t upperBound = 100.30;
  Int_t BinNumber=20;

  /* SIG specific variables */
  Double_t sigE=92.29;
  Double_t capturingRate=0.61; // capturing rate of Aluminum

  /* RMC specific vairables */
  Double_t PairCreationProb=0.97;
  Double_t VinTProb=0.0058; // Probability that vertex is in target
  Double_t ECut=0.016;      // Probability that positron has energy hihgher than 90.30

  ///////////////////////////
  // Set SIG and RMC Model //
  ///////////////////////////

  RooRealVar x("x","x",0,upperBound-lowerBound);

  // SIG component

  /* Landau distribuition determined by checking smeard energy from muon stopping target */
  Double_t sigNum=protonNum*muonStopRate*capturingRate*SIGBr*eff;

  RooRealVar landauMean("landauMean","mean of Landau",0.447);
  RooRealVar landauVar("landauVar","variance of Landau", 0.112);  
  RooRealVar sigEVal("sigEVal","signal energy", 92.29);
  RooRealVar lowerBoundVal("lowerBoundVal","lowerBoundVal", 90.30);

  //RooLandau sigPdf("landau", "landau p.d.f.",x, landauMean, landauVar);
  //RooFormularVar x_modi("x_flipping", "shift-x",RooArgSet(shift,x));

  //RooAbsPdf* sigPdf = bindPdf("sig",ROOT::Math::landau_pdf,x,landauMean,landauVar);

  RooGenericPdf sigPdf("sig","sig","Landau(-x+92.29-90.30,landauMean,landauVar)",RooArgSet(x,sigEVal, lowerBoundVal,landauMean,landauVar));  
  RooRealVar shift("shift","shift",sigE-lowerBound);
  RooFormulaVar x_shift("x_shift", "x+shift",RooArgSet(x,shift));  
  //RooGenericPdf sigPdf("sigPdf","signal p.d.f",x,x_shift

  //TF1 *landau = new TF1("landau","92.29-TMath::Landau(x,[0],[1])",0,100);
  //landau->SetParameters(0.447,0.112);
  //RooAbsReal *sigPdf=bindFunction(landau,x);


  // RMC component
  
  /* Nedd to determine exponential variable of rmc distribution */
  Double_t rmcNum=protonNum*muonStopRate*RMCBr*PairCreationProb*VinTProb*ECut*eff;
  RooRealVar expConst("expConst", "exponential component of RMC", -1, -1., 0.);
  //RooExponential rmcPdf_tmp("rmc", "rmc p.d.f.", x, expConst);
  RooRealVar a0("a0","a0",0);
  RooRealVar a1("a1","a1",0);
  RooRealVar a2("a2","a2",0);
  RooRealVar a3("a3","a3",0);
  RooRealVar a4("a4","a4",0);
  RooRealVar gausMean("gausMean", "mean of gaussian", -10,10);
  RooRealVar gausVar("gausVar", "variance of gaussian", 0, 0., 10.);

  //RooGenericPdf rmcPdf_tmp("rmc", "a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x", RooArgSet(x,a0,a1,a2,a3,a4));
  RooGaussian rmcPdf_tmp("rmc", "rmc gaussian", x, gausMean,gausVar);

  std::string rootPath="./";
  std::string rootFile="extrmc_1e7.root";

  TFile *f = TFile::Open(TString(rootPath+rootFile));
  TTree *t = (TTree*)f->Get("trdata");
  TH1F* rmc_hist = new TH1F("rmc_hist","rmc_hist", BinNumber, 0, upperBound-lowerBound);

  Float_t Pairep_genTrE;
  t->SetBranchAddress("Pairep_genTrE",&Pairep_genTrE);
  for (Int_t i_evt=0; i_evt<1000000; i_evt++){
    t->GetEntry(i_evt);
    if (Pairep_genTrE>lowerBound && Pairep_genTrE<upperBound){
      rmc_hist->Fill(Pairep_genTrE-lowerBound);
    }
  }

  RooRealVar x_rmc("x_rmc","x_rmc",0,upperBound-lowerBound);
  RooDataSet rmc_data("rmc","rmc", x_rmc, Import(*rmc_hist));
  rmcPdf_tmp.fitTo(rmc_data, Range(0,upperBound-lowerBound));

  //RooRealVar expConst_fit("expConst_fit", "expConst_fit",expConst.getValV());
  RooRealVar gausMean_fit("expConst_fit", "expConst_fit",gausMean.getValV());
  RooRealVar gausVar_fit("expConst_fit", "expConst_fit",gausVar.getValV());
  //RooExponential rmcPdf("rmc", "rmc p.d.f.", x, expConst_fit);
  RooGaussian rmcPdf("rmc", "rmc p.d.f.", x, gausMean_fit, gausVar_fit);
  

  std::cout << "Number of Signal: " << sigNum << "  Number of RMC: " << rmcNum << std::endl;

  /////////////////////////
  // Set Composite Model //
  /////////////////////////
  
  RooRealVar sigFrac("rmcFrac", "Fraction of RMC", Double_t(sigNum)/(sigNum+rmcNum));
  RooRealVar rmcFrac("rmcFrac", "Fraction of RMC", Double_t(rmcNum)/(sigNum+rmcNum));

  /* Binning */
  RooBinning bins(BinNumber,0,upperBound-lowerBound);

  /* Make Composite Model */
  RooAddPdf modelPdf("model", "model", RooArgList(sigPdf, rmcPdf),sigFrac);
  RooDataSet *modelData=modelPdf.generate(x, sigNum+rmcNum);
  
  RooDataSet *rmc_x = rmcPdf.generate(x,rmcNum);
  RooPlot *rmcFrame = x.frame(Title("RMC PDF"));
  rmc_x->plotOn(rmcFrame,Binning(bins));
  rmcPdf.plotOn(rmcFrame);

  //RooPlot *rmc_dataFrame = rmc_
  //rmc_data.plotOn(rmcFrame,Binning(bins));
  
  /* Plot Composite Model */
  
  RooPlot* modelFrame = x.frame(Title("Composite PDF"));
  modelData->plotOn(modelFrame,Binning(bins));
  modelPdf.plotOn(modelFrame, Components(rmcPdf), LineStyle(9), LineColor(kBlue));
  modelPdf.plotOn(modelFrame, LineStyle(kDashed), LineColor(kRed));
  modelPdf.plotOn(modelFrame, Components(sigPdf), LineStyle(9), LineColor(kBlack));
  
  TCanvas* c = new TCanvas("Composite_Model","Composite_Model",1200,400) ;            
  c->Divide(3);
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; modelFrame->GetYaxis()->SetTitleOffset(1.6) ; modelFrame->Draw() ;    
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; rmcFrame->GetYaxis()->SetTitleOffset(1.6) ; rmcFrame->Draw() ;
  c->cd(3) ; rmc_hist->Draw();
  

}
