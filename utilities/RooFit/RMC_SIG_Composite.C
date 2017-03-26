#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLandau.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooBinning.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooChi2Var.h"
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

  RooRealVar x("x","x",lowerBound,upperBound);

  // SIG component

  /* Landau distribuition determined by checking smeard energy from muon stopping target */
  Double_t sigNum=protonNum*muonStopRate*capturingRate*SIGBr*eff;

  RooRealVar landauMean("landauMean","mean of Landau",0.447);
  RooRealVar landauVar("landauVar","variance of Landau", 0.112);  
  RooRealVar sigEVal("sigEVal","signal energy", 92.29);
  RooRealVar lowerBoundVal("lowerBoundVal","lowerBoundVal", 90.30);

  RooGenericPdf sigPdf("sig","sig","Landau(-x+sigEVal,landauMean,landauVar)",RooArgSet(x,sigEVal, lowerBoundVal,landauMean,landauVar));  

  //TF1 *landau = new TF1("landau","92.29-TMath::Landau(x,[0],[1])",0,100);
  //landau->SetParameters(0.447,0.112);
  //RooAbsReal *sigPdf=bindFunction(landau,x);

  // RMC component
  
  /* Nedd to determine exponential variable of rmc distribution */
  Double_t rmcNum=protonNum*muonStopRate*RMCBr*PairCreationProb*VinTProb*ECut*eff;
  RooRealVar expConst("expConst", "exponential component of RMC", -1, -1., 0.);
  RooRealVar a0("a0","a0",0);
  RooRealVar a1("a1","a1",0);
  RooRealVar a2("a2","a2",0);
  RooRealVar a3("a3","a3",0);
  RooRealVar a4("a4","a4",0);
  RooRealVar gausMean("gausMean", "mean of gaussian", -30,-50.,-30.);
  RooRealVar gausVar("gausVar", "variance of gaussian", 10, 5., 15.);

  //RooGenericPdf rmcPdf_tmp("rmc", "a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x", RooArgSet(x,a0,a1,a2,a3,a4));
  RooRealVar x_fit("x_fit","x_fit",0,upperBound-lowerBound);
  RooGaussian rmcPdf_tmp("rmc", "rmc gaussian", x_fit, gausMean,gausVar);

  std::string rootPath="./";
  std::string rootFile="extrmc_1e7.root";

  TFile *f = TFile::Open(TString(rootPath+rootFile));
  TTree *t = (TTree*)f->Get("trdata");
  TH1F* rmc_hist = new TH1F("rmc_hist","rmc_hist", BinNumber, 0, upperBound-lowerBound);

  Float_t Pairep_genTrE;
  t->SetBranchAddress("Pairep_genTrE",&Pairep_genTrE);
  for (Int_t i_evt=0; i_evt<100000; i_evt++){
    t->GetEntry(i_evt);
    if (Pairep_genTrE>lowerBound && Pairep_genTrE<upperBound){
      rmc_hist->Fill(Pairep_genTrE-lowerBound);
    }
  }

  RooRealVar x_rmc("x_rmc","x_rmc",0,upperBound-lowerBound);
  RooDataSet rmc_data("rmc","rmc", x_rmc, Import(*rmc_hist));
  rmcPdf_tmp.fitTo(rmc_data);

  RooRealVar gausMeanFit("gausMeanFit", "gausMeanFit",gausMean.getValV()+lowerBound);
  RooRealVar gausVarFit("gausVarFit", "gausVarFit",gausVar.getValV());
  RooGaussian rmcPdf("rmc", "rmc p.d.f.", x, gausMeanFit, gausVarFit);
  
  std::cout << "RMC Gaus Mean: " <<gausMean.getValV() << "  RMC Gaus Variance: " << gausVarFit.getValV() << std::endl;
  std::cout << "Number of Signal: " << sigNum << "  Number of RMC: " << rmcNum << std::endl;


  /////////////////////////
  // Set Composite Model //
  /////////////////////////
  
  RooRealVar sigFrac("rmcFrac", "Fraction of RMC", Double_t(sigNum)/(sigNum+rmcNum));
  RooRealVar rmcFrac("rmcFrac", "Fraction of RMC", Double_t(rmcNum)/(sigNum+rmcNum));

  /* Binning */
  RooBinning bins(BinNumber,lowerBound,upperBound);

  /* Make Composite Model */
  RooAddPdf modelPdf("model", "model", RooArgList(sigPdf, rmcPdf),sigFrac);
  RooDataSet *modelData=modelPdf.generate(x, (sigNum+rmcNum)/1000);
   
  /* Plot Composite Model */
  
  RooPlot* modelFrame = x.frame(Title("Composite PDF"));
  modelData->plotOn(modelFrame,Binning(bins));
  modelPdf.plotOn(modelFrame, Components(rmcPdf), LineStyle(9), LineColor(kBlue));
  modelPdf.plotOn(modelFrame, LineStyle(kDashed), LineColor(kRed));
  //modelPdf.plotOn(modelFrame, Components(sigPdf), LineStyle(9), LineColor(kBlack));

  TCanvas* c = new TCanvas("Composite_Model","Composite_Model",800,400);  
  c->Divide(2);
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; modelFrame->GetYaxis()->SetTitleOffset(1.6) ; modelFrame->Draw() ;    
  //c->cd(2) ; gPad->SetLeftMargin(0.15) ; rmcFrame->GetYaxis()->SetTitleOffset(1.6) ; rmcFrame->Draw() ;
  c->cd(2) ; rmc_hist->Draw();
  

}
