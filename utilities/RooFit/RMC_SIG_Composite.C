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
#include "RooAbsArg.h"
#include "RooFitResult.h"
#include "RooFit.h"
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
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"

void DoHypothesisTest(RooWorkspace*);
void MakePlots(RooWorkspace*);


using namespace RooFit;
using namespace RooStats;

void RMC_SIG_Composite(){
  
  ///////////////////
  // Set Variables //
  ///////////////////

  /* branching ratios */
  Double_t SIGBr=8*TMath::Power(10,-13); // Set to current limit on mu->e+ conversion
  Double_t RMCBr=6.22*TMath::Power(10,-7);

  /* common variables */
  Double_t eff=0.01;
  Double_t protonNum=TMath::Power(10,21);
  Double_t muonStopRate=TMath::Power(10,-3);
  Double_t lowerBound = 90.30;
  Double_t upperBound = 95.80;
  Int_t BinNumber=11;

  /* SIG specific variables */
  Double_t sigE=92.29;
  Double_t capturingRate=0.61; // capturing rate of Aluminum

  /* RMC specific vairables */
  Double_t PairCreationProb=0.96;
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
  RooRealVar sigEVal("sigEVal","signal energy", sigE);
  RooRealVar lowerBoundVal("lowerBoundVal","lowerBoundVal", lowerBound);

  RooGenericPdf sigPdf("sig","sig","Landau(-x+sigEVal,landauMean,landauVar)",RooArgSet(x,sigEVal, lowerBoundVal,landauMean,landauVar));  

  //TF1 *landau = new TF1("landau","92.29-TMath::Landau(x,[0],[1])",0,100);
  //landau->SetParameters(0.447,0.112);
  //RooAbsReal *sigPdf=bindFunction(landau,x);

  // RMC component
  
  /* Nedd to determine exponential variable of rmc distribution */
  Double_t rmcNum=protonNum*muonStopRate*RMCBr*PairCreationProb*VinTProb*ECut*eff; //*2; //factor 2 from internal conversion
  RooRealVar expConst("expConst", "exponential component of RMC", -0.37, -0.39, -0.35);
  //RooRealVar expConst("expConst", "exponential component of RMC", -0.45, -0.60, -0.34);
  RooRealVar linConst("linConst", "linear compenent of RMC", 0.1, 0.1 , 0.2);

  //RooRealVar a0("a0","a0",0);
  //RooRealVar a1("a1","a1",0);
  //RooRealVar a2("a2","a2",0);
  //RooRealVar a3("a3","a3",0);
  //RooRealVar a4("a4","a4",0);
  //RooRealVar gausMean("gausMean", "mean of gaussian", -2,-5.,-1.);
  //RooRealVar gausVar("gausVar", "variance of gaussian", 2.9999, 0., 5.);

  RooRealVar x_fit("x_fit","x_fit",0,upperBound-lowerBound);
  RooExponential rmcPdf_tmp("rmc", "rmc exponential", x_fit, expConst);
  //RooGenericPdf rmcPdf_tmp("rmc", "exp(x_fit*expConst)/(x_fit+linConst)", RooArgSet(x_fit,expConst,linConst));

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

  RooGenericPdf rmcPdf("rmc","rmc","exp((x-lowerBoundVal)*expConst)",RooArgSet(x,lowerBoundVal,expConst));    
//RooGenericPdf rmcPdf("rmc","rmc","exp((x-lowerBoundVal)*expConst)/((x-lowerBoundVal)+linConst)",RooArgSet(x,lowerBoundVal,expConst,linConst));

  //std::cout << "RMC Gaus Mean: " <<gausMean.getValV() << "  RMC Gaus Variance: " << gausVarFit.getValV() << std::endl;

  /////////////////////////
  // Set Composite Model //
  /////////////////////////
  
  RooRealVar sigFrac("sigFrac", "Fraction of RMC", Double_t(sigNum)/(sigNum+rmcNum),0.,1.);
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


  RooAbsReal* rmcInt = rmcPdf.createIntegral(x,NormSet(x),Range("ECut"));
  RooAbsReal* sigInt = sigPdf.createIntegral(x,NormSet(x),Range("ECut"));
  std::cout << "-------------------------------------------------" << std::endl;
  expConst.Print();
  std::cout << "Number of Signal: " << sigNum << "  Number of RMC: " << rmcNum << std::endl;
  x.setRange("ECut",91.7,92.0);
  std::cout << "Number of Signal in ECut Boundary " << sigInt->getValV()*sigNum << std::endl;
  std::cout << "Number of RMC    in ECut Boundary " << rmcInt->getValV()*rmcNum << std::endl;
  std::cout << "-------------------------------------------------\n\n" << std::endl;


  /////////////////////////////////
  // Draw MC Data & Fit Function //
  /////////////////////////////////

  
  TCanvas* c2 = new TCanvas("Composite_Model_witHistogram","Composite Model with Histogram",600,600); 
  c2->cd() ;

  RooPlot* modelFrame2 = x.frame(Title("Composite PDF"));  
  RooDataSet *modelData2=modelPdf.generate(x, rmc_hist2->GetEntries());
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
  rmc_graph->GetYaxis()->SetRangeUser(290,3060);
  rmc_graph->GetYaxis()->SetAxisColor(kRed);
  rmc_graph->Draw();
  pad1_2->Draw();
  pad1_2->cd();  
  rmc_hist2->Draw();


  ///////////////////////////////
  // RooStats Likelihood Method //
  ///////////////////////////////

  /* Set up Model with ProfileLikelihoodCalculator */
    
   
  ModelConfig model;
  RooWorkspace* wks = new RooWorkspace("wks");
  wks->import(modelPdf);
  
  wks->import(*modelData, Rename("data"));
  model.SetWorkspace(*wks);
  model.SetPdf("model");
  
  
  ProfileLikelihoodCalculator plc;
  plc.SetData(*(wks->data("data")));
  RooRealVar* fsig = wks->var("sigFrac");
  RooArgSet poi(*fsig);
  RooArgSet *nullParams = (RooArgSet*)poi.snapshot();
  nullParams->setRealValue("sigFrac",0);

  plc.SetModel(model);
  plc.SetNullParameters(*nullParams);
  //plc.SetNullParameters(*nullParams);                                            
  plc.SetModel(model);
  // NOTE: using snapshot will import nullparams                                   
  // in the WS and merge with existing "mu"                                        
  // model.SetSnapshot(*nullParams);                                               
  
  //use instead setNuisanceParameters                                              
  plc.SetNullParameters( *nullParams);
  
  // We get a HypoTestResult out of the calculator, and we can query it.           
  HypoTestResult* htr = plc.GetHypoTest();
  cout << "-------------------------------------------------" << endl;
  cout << "The p-value for the null is " << htr->NullPValue() << endl;
  cout << "Corresponding to a significance of " << htr->Significance() << endl;
  cout << "-------------------------------------------------\n\n" << endl;
    
  //MakePlots(wks);
  delete wks;  

  
}
