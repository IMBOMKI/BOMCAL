#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <tuple>


void AnalyzeOutput(){
  Int_t NumOfEvt=100000;
  Int_t bw_ini=3;
  Int_t bw_fin=10;
  Double_t bw[bw_fin-bw_ini+1];
  Double_t Acceptance_em[bw_fin-bw_ini+1];
  Double_t Acceptance_ep[bw_fin-bw_ini+1];
  Double_t MisCIDRate_em[bw_fin-bw_ini+1];
  Double_t MisCIDRate_ep[bw_fin-bw_ini+1];

  Double_t RecoRate_em[bw_fin-bw_ini+1];
  Double_t RecoRate_ep[bw_fin-bw_ini+1];
  Double_t RecoRate_Accep_em[bw_fin-bw_ini+1];
  Double_t RecoRate_Accep_ep[bw_fin-bw_ini+1];
  memset (Acceptance_em, 0, sizeof (Acceptance_em));
  memset (Acceptance_ep, 0, sizeof (Acceptance_ep));
  memset (RecoRate_em, 0, sizeof (RecoRate_em));
  memset (RecoRate_ep, 0, sizeof (RecoRate_ep));

  std::string dir = "./FindedEvents_niter3_nPt300/"; 
  ///////////////////////////////////////////////////////////////////////////////
  //
  //
  // e- Part
  //
  //
  ///////////////////////////////////////////////////////////////////////////////
  
  for (Int_t i_bw=bw_ini; i_bw<bw_fin+1; i_bw++){
    
    bw[i_bw-bw_ini] = i_bw;
    std::string fileName_em = "bw"+std::to_string(i_bw)+"_find_trig_em104_onlyPrimary.root";
    TFile *f_em = TFile::Open(TString(dir+fileName_em));
    TTree *t_em = (TTree*)f_em->Get("trdata");

    Int_t nCALCDCHit;
    Int_t nRecoHit;
    Int_t RecoMaxWireLayerId;
    Bool_t Reco_ifCL3;
    Bool_t ifSingleTurn;
    Bool_t ifMultiTurn;
    Double_t fittedR_even;
    Double_t fittedR_odd;
    Double_t drEvenToOdd;
    Int_t RecoCharge;

    t_em->SetBranchAddress("nCALCDCHit", &nCALCDCHit);    
t_em->SetBranchAddress("nRecoHit", &nRecoHit);
    t_em->SetBranchAddress("RecoMaxWireLayerId", &RecoMaxWireLayerId);
    t_em->SetBranchAddress("Reco_ifCL3", &Reco_ifCL3);
    t_em->SetBranchAddress("ifSingleTurn", &ifSingleTurn);
    t_em->SetBranchAddress("ifMultiTurn", &ifMultiTurn);
    t_em->SetBranchAddress("fittedR_even", &fittedR_even);
    t_em->SetBranchAddress("fittedR_odd", &fittedR_odd);
    t_em->SetBranchAddress("drEvenToOdd", &drEvenToOdd);
    t_em->SetBranchAddress("RecoCharge", &RecoCharge);

    Double_t TotalHit=0;
    Double_t Total_RecoHit=0;
    Double_t TotalHit_Accep=0;
    Double_t Total_RecoHit_Accep=0;

    for (Int_t i_evt=0; i_evt<t_em->GetEntries(); i_evt++){      
      t_em->GetEntry(i_evt);
      
      // Acceptance && Charge MisIdentification

      if (Reco_ifCL3==1 && RecoMaxWireLayerId>=4 && nRecoHit>=30 && drEvenToOdd<14 && fittedR_even<40 && fittedR_odd<40){
	Acceptance_em[i_bw-bw_ini]++;
	TotalHit_Accep+=nCALCDCHit;
	Total_RecoHit_Accep+=nRecoHit;

	if (RecoCharge==1){
	  MisCIDRate_em[i_bw-bw_ini]++;
	}
      }

      // Reco Rate

      TotalHit+=nCALCDCHit;
      Total_RecoHit+=nRecoHit;
    }

    MisCIDRate_em[i_bw-bw_ini]=MisCIDRate_em[i_bw-bw_ini]/Double_t(Acceptance_em[i_bw-bw_ini])*100;
    Acceptance_em[i_bw-bw_ini]=Acceptance_em[i_bw-bw_ini]/NumOfEvt*100;
    RecoRate_em[i_bw-bw_ini]=Double_t(Total_RecoHit)/(TotalHit)*100.0;
    RecoRate_Accep_em[i_bw-bw_ini]=Double_t(Total_RecoHit_Accep)/(TotalHit_Accep)*100.0;
    
    f_em->cd();
    f_em->Close();
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  //
  //
  // e+ Part
  //
  //
  ///////////////////////////////////////////////////////////////////////////////
  
  for (Int_t i_bw=bw_ini; i_bw<bw_fin+1; i_bw++){
    
    bw[i_bw-bw_ini] = i_bw;
    std::string fileName_ep = "bw"+std::to_string(i_bw)+"_find_trig_ep93_onlyPrimary.root";
    TFile *f_ep = TFile::Open(TString(dir+fileName_ep));
    TTree *t_ep = (TTree*)f_ep->Get("trdata");

    Int_t nCALCDCHit;
    Int_t nRecoHit;
    Int_t RecoMaxWireLayerId;
    Bool_t Reco_ifCL3;
    Bool_t ifSingleTurn;
    Bool_t ifMultiTurn;
    Double_t fittedR_even;
    Double_t fittedR_odd;
    Double_t drEvenToOdd;
    Int_t RecoCharge;

    t_ep->SetBranchAddress("nCALCDCHit", &nCALCDCHit);
    t_ep->SetBranchAddress("nRecoHit", &nRecoHit);
    t_ep->SetBranchAddress("RecoMaxWireLayerId", &RecoMaxWireLayerId);
    t_ep->SetBranchAddress("Reco_ifCL3", &Reco_ifCL3);
    t_ep->SetBranchAddress("ifSingleTurn", &ifSingleTurn);
    t_ep->SetBranchAddress("ifMultiTurn", &ifMultiTurn);
    t_ep->SetBranchAddress("fittedR_even", &fittedR_even);
    t_ep->SetBranchAddress("fittedR_odd", &fittedR_odd);
    t_ep->SetBranchAddress("drEvenToOdd", &drEvenToOdd);
    t_ep->SetBranchAddress("RecoCharge", &RecoCharge);

    Double_t TotalHit=0;
    Double_t Total_RecoHit=0;
    Double_t TotalHit_Accep=0;
    Double_t Total_RecoHit_Accep=0;


    for (Int_t i_evt=0; i_evt<t_ep->GetEntries(); i_evt++){      
      t_ep->GetEntry(i_evt);
      
      // Acceptance

      if (Reco_ifCL3==1 && RecoMaxWireLayerId>=4 && nRecoHit>=30 && drEvenToOdd<14 && fittedR_even<45 && fittedR_odd<45){
	Acceptance_ep[i_bw-bw_ini]++;
	TotalHit_Accep+=nCALCDCHit;
	Total_RecoHit_Accep+=nRecoHit;

	if (RecoCharge==-1){
	  MisCIDRate_ep[i_bw-bw_ini]++;
	}
      }

      // Reco Rate

      TotalHit+=nCALCDCHit;
      Total_RecoHit+=nRecoHit;
    }

    //std::cout <<     MisCIDRate_ep[i_bw-bw_ini] << std::endl;

    MisCIDRate_ep[i_bw-bw_ini]=MisCIDRate_ep[i_bw-bw_ini]/Double_t(Acceptance_ep[i_bw-bw_ini])*100;
    Acceptance_ep[i_bw-bw_ini]=Acceptance_ep[i_bw-bw_ini]/NumOfEvt*100;
    RecoRate_ep[i_bw-bw_ini]=Double_t(Total_RecoHit)/(TotalHit)*100.0;
    RecoRate_Accep_ep[i_bw-bw_ini]=Double_t(Total_RecoHit_Accep)/(TotalHit_Accep)*100.0;

    std::cout <<  Acceptance_em[i_bw-bw_ini] << "   " << MisCIDRate_em[i_bw-bw_ini] << "   " << RecoRate_em[i_bw-bw_ini]  << "   " <<  Acceptance_ep[i_bw-bw_ini] << "   " << MisCIDRate_ep[i_bw-bw_ini] << "   " << RecoRate_ep[i_bw-bw_ini] << std::endl;


    f_ep->cd();
    f_ep->Close();
  }
  

  TCanvas *c_div = new TCanvas("canvas", "canvas", 2000, 500);
  c_div->Divide(4,1);
  TPad *pad1_1 = new TPad("pad1","",0,0,1,1);
  TPad *pad1_2 = new TPad("pad2","",0,0,1,1);
  pad1_1->SetFillStyle(4000);
  pad1_1->SetFrameFillStyle(0);
  pad1_2->SetFillStyle(4000);
  pad1_2->SetFrameFillStyle(0);

  TPad *pad2_1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2_2 = new TPad("pad2","",0,0,1,1);
  pad2_1->SetFillStyle(4000);
  pad2_1->SetFrameFillStyle(0);
  pad2_2->SetFillStyle(4000);
  pad2_2->SetFrameFillStyle(0);

  TPad *pad3_1 = new TPad("pad1","",0,0,1,1);  // Reco Rate for All Events
  TPad *pad3_2 = new TPad("pad2","",0,0,1,1);
  pad3_2->SetFillStyle(4000);
  pad3_2->SetFrameFillStyle(0);

  TPad *pad4_1 = new TPad("pad1","",0,0,1,1);  // Reco Rate for only Accepted Events
  TPad *pad4_2 = new TPad("pad2","",0,0,1,1);
  pad4_2->SetFillStyle(4000);
  pad4_2->SetFrameFillStyle(0);

  TMultiGraph *mg1 = new TMultiGraph(); 
  TGraph *grAcceptance_em= new TGraph(bw_fin-bw_ini+1, bw, Acceptance_em);
  TGraph *grMisCIDRate_em= new TGraph(bw_fin-bw_ini+1, bw, MisCIDRate_em);
  TGraph *grRecoRate_em= new TGraph(bw_fin-bw_ini+1, bw, RecoRate_em);
  TGraph *grRecoRate_Accep_em= new TGraph(bw_fin-bw_ini+1, bw, RecoRate_Accep_em);

  TGraph *grAcceptance_ep= new TGraph(bw_fin-bw_ini+1, bw, Acceptance_ep);
  TGraph *grMisCIDRate_ep= new TGraph(bw_fin-bw_ini+1, bw, MisCIDRate_ep);
  TGraph *grRecoRate_ep= new TGraph(bw_fin-bw_ini+1, bw, RecoRate_ep);
  TGraph *grRecoRate_Accep_ep= new TGraph(bw_fin-bw_ini+1, bw, RecoRate_Accep_ep);

  ///////////////////////////////////////////////////////////////////////
  c_div->cd(1);
  pad1_1->Draw();
  pad1_1->cd();
  grAcceptance_em->SetTitle("Acceptance");
  grAcceptance_em->SetMarkerStyle(20);
  grAcceptance_em->SetMarkerColor(4);
  grAcceptance_em->GetYaxis()->SetRangeUser(9,20);
  grAcceptance_em->GetYaxis()->SetAxisColor(4);
  grAcceptance_em->GetYaxis()->SetTitle("Acceptance of e- (%)");
  grAcceptance_em->GetXaxis()->SetTitle("Band Width of Circle");
  grAcceptance_em->Draw("ALP");
  
  pad1_2->Draw();
  pad1_2->cd();
  grAcceptance_ep->SetMarkerStyle(20);
  grAcceptance_ep->SetMarkerColor(2);
  grAcceptance_ep->GetYaxis()->SetRangeUser(2,3.5);
  grAcceptance_ep->GetYaxis()->SetAxisColor(2);
  grAcceptance_ep->GetYaxis()->SetTitle("Acceptance of e+ (%)");
  grAcceptance_ep->Draw("LPAY+");

  /////////////////////////////////////////////////////////////////////
  
  c_div->cd(2);
  pad2_1->Draw();
  pad2_1->cd();
  grMisCIDRate_em->SetTitle("Mis Charge Identification Rate");
  grMisCIDRate_em->SetMarkerStyle(20);
  grMisCIDRate_em->SetMarkerColor(4);
  grMisCIDRate_em->GetYaxis()->SetRangeUser(0,1);
  grMisCIDRate_em->GetYaxis()->SetAxisColor(4);
  grMisCIDRate_em->GetYaxis()->SetTitle("Mis CID of e- as e+ (%)");
  grMisCIDRate_em->GetXaxis()->SetTitle("Band Width of Circle");
  grMisCIDRate_em->Draw("ALP");
  
  pad2_2->Draw();
  pad2_2->cd();
  grMisCIDRate_ep->SetMarkerStyle(20);
  grMisCIDRate_ep->SetMarkerColor(2);
  grMisCIDRate_ep->GetYaxis()->SetRangeUser(5,15);
  grMisCIDRate_ep->GetYaxis()->SetAxisColor(2);
  grMisCIDRate_ep->GetYaxis()->SetTitle("Mis CID of e+ as e- (%)");
  grMisCIDRate_ep->Draw("LPAY+");
  
  ////////////////////////////////////////////////////////////////////

  c_div->cd(3);
  pad3_1->Draw();
  pad3_1->cd();
  grRecoRate_em->SetTitle("Signal Hits Recognition Rate");
  grRecoRate_em->SetMarkerStyle(20);
  grRecoRate_em->SetMarkerColor(38);
  grRecoRate_em->GetYaxis()->SetRangeUser(93,100);
  grRecoRate_em->GetYaxis()->SetTitle("Recognition Rate (%)");
  grRecoRate_em->GetXaxis()->SetTitle("Band Width of Circle");
  grRecoRate_em->Draw("ALP");
  
  pad3_2->Draw();
  pad3_2->cd();
  grRecoRate_ep->SetTitle("Signal Hits Recognition Rate");
  grRecoRate_ep->SetMarkerStyle(20);
  grRecoRate_ep->SetMarkerColor(46);
  grRecoRate_ep->GetYaxis()->SetRangeUser(93,100);
  grRecoRate_ep->Draw("LPA");

  //////////////////////////////////

  c_div->cd(4);
  pad4_1->Draw();
  pad4_1->cd();
  grRecoRate_Accep_em->SetTitle("Signal Hits Recognition Rate for Accepted Events");
  grRecoRate_Accep_em->SetMarkerStyle(20);
  grRecoRate_Accep_em->SetMarkerColor(38); // electron(e-) Marker is light blue
  grRecoRate_Accep_em->GetYaxis()->SetRangeUser(93,100);
  grRecoRate_Accep_em->GetYaxis()->SetTitle("Recognition Rate (%)");
  grRecoRate_Accep_em->GetXaxis()->SetTitle("Band Width of Circle");
  grRecoRate_Accep_em->Draw("ALP");
  
  pad4_2->Draw();
  pad4_2->cd();
  grRecoRate_Accep_ep->SetTitle("Signal Hits Recognition Rate for Accepted Events");
  grRecoRate_Accep_ep->SetMarkerStyle(20);  
  grRecoRate_Accep_ep->SetMarkerColor(46); // positron(e+) Marker is light red
  grRecoRate_Accep_ep->GetYaxis()->SetRangeUser(93,100);
  grRecoRate_Accep_ep->Draw("LPA");



}

