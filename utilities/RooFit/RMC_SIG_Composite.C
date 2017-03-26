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
#include "TGraph.h"
#include "TPad.h"
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
  RooRealVar expConst("expConst", "exponential component of RMC", -0.40, -0.50, -0.34);
  //RooRealVar expConst("expConst", "exponential component of RMC", -5, -5., -4.);
  RooRealVar a0("a0","a0",0);
  RooRealVar a1("a1","a1",0);
  RooRealVar a2("a2","a2",0);
  RooRealVar a3("a3","a3",0);
  RooRealVar a4("a4","a4",0);
  RooRealVar gausMean("gausMean", "mean of gaussian", -2,-5.,-1.);
  RooRealVar gausVar("gausVar", "variance of gaussian", 2.9999, 0., 5.);

  //RooGenericPdf rmcPdf_tmp("rmc", "a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x", RooArgSet(x,a0,a1,a2,a3,a4));
  RooRealVar x_fit("x_fit","x_fit",0,upperBound-lowerBound);
  //RooGaussian rmcPdf_tmp("rmc", "rmc gaussian", x_fit, gausMean,gausVar);
  RooExponential rmcPdf_tmp("rmc", "rmc exponential", x_fit, expConst);

  std::string rootPath="./";
  std::string rootFile="extrmc_1e7_second.root";

  TFile *f = TFile::Open(TString(rootPath+rootFile));
  TTree *t = (TTree*)f->Get("trdata");
  TH1F* rmc_hist = new TH1F("rmc_hist","rmc_hist", BinNumber, 0, upperBound-lowerBound);

  Float_t Pairep_genTrE;
  t->SetBranchAddress("Pairep_genTrE",&Pairep_genTrE);

  TH1F* rmc_hist2 = new TH1F("rmc_hist2","rmc_hist2", BinNumber, lowerBound, upperBound);
  for (Int_t i_evt=0; i_evt<1000000; i_evt++){
    t->GetEntry(i_evt);
    if (Pairep_genTrE>lowerBound && Pairep_genTrE<upperBound){
      rmc_hist->Fill(Pairep_genTrE-lowerBound);
      rmc_hist2->Fill(Pairep_genTrE);
    }
  }

  RooRealVar x_rmc("x_rmc","x_rmc",0,upperBound-lowerBound);
  RooDataSet rmc_data("rmc","rmc", x_rmc, Import(*rmc_hist));
  rmcPdf_tmp.fitTo(rmc_data, Range(0,upperBound-lowerBound));

  //RooRealVar gausMeanFit("gausMeanFit", "gausMeanFit",gausMean.getValV()+lowerBound);
  //RooRealVar gausVarFit("gausVarFit", "gausVarFit",gausVar.getValV());
  //RooGaussian rmcPdf("rmc", "rmc p.d.f.", x, gausMeanFit, gausVarFit);

  //RooRealVar expConstFit("expConstFit", "expConstFit",expConst.getValV());
  //RooExponential rmcPdf("rmc", "rmc p.d.f.", x, expConst);
  RooGenericPdf rmcPdf("rmc","rmc","exp((x-lowerBoundVal)*expConst)",RooArgSet(x,lowerBoundVal,expConst));    

  //std::cout << "RMC Gaus Mean: " <<gausMean.getValV() << "  RMC Gaus Variance: " << gausVarFit.getValV() << std::endl;
  expConst.Print();
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
  RooDataSet *modelData=modelPdf.generate(x, (sigNum+rmcNum));
   
  /* Plot Composite Model */
  RooPlot* modelFrame = x.frame(Title(""));
  modelFrame->SetTitle("");
  modelData->plotOn(modelFrame,Binning(bins));
  modelPdf.plotOn(modelFrame, Components(rmcPdf), LineStyle(9), LineColor(kBlue));
  modelPdf.plotOn(modelFrame, LineStyle(kDashed), LineColor(kRed));
  //modelPdf.plotOn(modelFrame, Components(sigPdf), LineStyle(9), LineColor(kBlack));

  TCanvas* c1 = new TCanvas("Composite_Model","Composite_Model",600,600);  
  c1->cd() ; 
  gPad->SetLeftMargin(0.15) ; 
  modelFrame->GetYaxis()->SetTitleOffset(1.8) ; 
  modelFrame->GetXaxis()->SetTitle("E_{e^{+}} (MeV)");
  modelFrame->Draw() ;    

  TCanvas* c2 = new TCanvas("Composite_Model_witHistogram","Composite Model with Histogram",600,600); 
  c2->cd() ;

  RooPlot* modelFrame2 = x.frame(Title("Composite PDF"));  
  RooDataSet *modelData2=modelPdf.generate(x, 15787);
  modelData2->plotOn(modelFrame2, Components(rmcPdf), LineStyle(9), LineColor(kBlue),Name("rmcGr"),Binning(bins));
  TGraph* rmc_graph = (TGraph*)modelFrame2->getObject( modelFrame2->numItems() - 1  );

  TPad *pad1_1 = new TPad("pad1","",0,0,1,1);
  TPad *pad1_2 = new TPad("pad2","",0,0,1,1);
  pad1_1->SetFillStyle(4000);
  pad1_1->SetFrameFillStyle(0);
  pad1_2->SetFillStyle(4000);
  pad1_2->SetFrameFillStyle(0);

  pad1_1->Draw();
  pad1_1->cd();
  rmc_graph->GetXaxis()->SetRangeUser(lowerBound,upperBound);
  rmc_graph->GetYaxis()->SetRangeUser(0,3070);
  rmc_graph->Draw();
  pad1_2->Draw();
  pad1_2->cd();  
  rmc_hist2->Draw();
}
