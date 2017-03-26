#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TChain.h"
#include "TGraph.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


void AnalyzeBunch(){

  TFile *f =TFile::Open("../2500BunchTrains.root");
  TTree *t = (TTree*)f->Get("trdata");

  Int_t NumOfTrack;
  Bool_t Tr_ifHitCDC[50000];
  Bool_t Tr_ifHitCRK[50000];
  Bool_t Tr_ifHitSTL[50000];
  Int_t PDGNum[50000];
  Float_t genTrE[50000];
  Float_t genTrPx[50000];
  Float_t genTrPy[50000];
  Float_t genTrPz[50000];
  Float_t genTrX[50000];
  Float_t genTrY[50000];
  Float_t genTrZ[50000];
  Float_t genTrT[50000];
  Int_t TrackId[50000];
  Int_t ParentId[50000];
  Bool_t Tr_ifStartedAtTarget[50000];
  Bool_t Tr_ifStoppedAtTarget[50000];

  t->SetBranchAddress("NumOfTrack",&NumOfTrack);
  t->SetBranchAddress("Tr_ifHitCDC",Tr_ifHitCDC);
  t->SetBranchAddress("Tr_ifHitCRK",Tr_ifHitCRK);
  t->SetBranchAddress("Tr_ifHitSTL",Tr_ifHitSTL);
  t->SetBranchAddress("PDGNum",PDGNum);
  t->SetBranchAddress("genTrE",genTrE);
  t->SetBranchAddress("genTrPx",genTrPx);
  t->SetBranchAddress("genTrPy",genTrPy);
  t->SetBranchAddress("genTrPz",genTrPz);
  t->SetBranchAddress("genTrX",genTrX);
  t->SetBranchAddress("genTrY",genTrY);
  t->SetBranchAddress("genTrZ",genTrZ);
  t->SetBranchAddress("genTrT",genTrT);
  t->SetBranchAddress("TrackId", TrackId);
  t->SetBranchAddress("ParentId", ParentId);
  t->SetBranchAddress("Tr_ifStartedAtTarget", Tr_ifStartedAtTarget);
  t->SetBranchAddress("Tr_ifStoppedAtTarget", Tr_ifStoppedAtTarget);

  Int_t NumOfBunch=2500;
  TString prim = "primary";
  TString em = "e-";
  TString ep = "e+";
  TString mu = "mu-";
  TString antimu = "mu+";
  TString photon = "photon";
  TString neutron = "neutron";
  TString pi0 = "pi0";
  TString others = "others";

  TH1F *h_ep = new TH1F("ep_hist","Energy of e+ passing detectors (CDC+STL+CRK)",110,0,110);  
  TH1I *h_ep_ori = new TH1I("origin_of_ep","Orgin of high energy e+ (>80 MeV)",10,0,10);
  TH1I *h_gamma_ori = new TH1I("origin_of_gamma","Orgin of high energy photon (>80 MeV)",10,0,10);
  TH1F *h_gamma_ene = new TH1F("energy_of_gamma","Photon energy (>80 MeV) coming from target",30,80,110);

  for (Int_t i_evt=0; i_evt<NumOfBunch; i_evt++){
    t->GetEntry(i_evt);
    std::cout << i_evt << "-th train is being processed..." << std::endl;  

    for (Int_t i_tr=0; i_tr<NumOfTrack; i_tr++){

      if(PDGNum[i_tr]==-11){

	// Count positrons passing through detectors (CDC+trigger Hodoscope)
	// and make energy histogram

	if(Tr_ifHitCDC[i_tr]==1 && Tr_ifHitCRK[i_tr]==1 && Tr_ifHitSTL[i_tr]==1){	  
	  h_ep->Fill(genTrE[i_tr]);
	}

	// Origin of high energy positron (Birth Energy > 80 MeV)

	if(genTrE[i_tr]>80){
	  if (ParentId[i_tr]==0){  // Primary
	      h_ep_ori->Fill(prim,1);
	  }
	  if (ParentId[i_tr]!=0){
	    for (Int_t i_tr2=0; i_tr2<NumOfTrack; i_tr2++){
	      if (TrackId[i_tr2]==ParentId[i_tr]){
		
		if (PDGNum[i_tr2]==11)           h_ep_ori->Fill(em,1);
		else if (PDGNum[i_tr2]==-11)     h_ep_ori->Fill(ep,1);
		else if (PDGNum[i_tr2]==13)      h_ep_ori->Fill(mu,1);
		else if (PDGNum[i_tr2]==-13)     h_ep_ori->Fill(antimu,1);
		else if (PDGNum[i_tr2]==22)      h_ep_ori->Fill(photon,1);
		else if (PDGNum[i_tr2]==2112)    h_ep_ori->Fill(neutron,1);
		//else if (PDGNum[i_tr2]==111)   h_ep_ori->Fill(pi0,1);
		else  {                    
                                  		 h_ep_ori->Fill(others,1); // others
		  std::cout << "Positron origin others: " << PDGNum[i_tr2] << std::endl;

		}
	      }
	    }
	  }
	}
      }

      
      // Analyze the origin of high energy photon independently
      
      if(PDGNum[i_tr]==22 && genTrE[i_tr]>80){
	
	  if (ParentId[i_tr]==0){  // Primary
	    h_gamma_ori->Fill(prim,1);
	  }
	  if (ParentId[i_tr]!=0){
	    for (Int_t i_tr2=0; i_tr2<NumOfTrack; i_tr2++){
	      if (TrackId[i_tr2]==ParentId[i_tr]){
		
		if (PDGNum[i_tr2]==11)          h_gamma_ori->Fill(em,1);
		else if (PDGNum[i_tr2]==-11)    h_gamma_ori->Fill(ep,1);   
		else if (PDGNum[i_tr2]==13)     h_gamma_ori->Fill(mu,1);   
		else if (PDGNum[i_tr2]==111)    h_gamma_ori->Fill(pi0,1);   
		else if   {                     
		                                h_gamma_ori->Fill(others,1);  // other
		  std::cout << "Gamma origin others: " << PDGNum[i_tr2] << std::endl;
		}
	      }
	    }
	  }
      }

      if(PDGNum[i_tr]==22 && genTrE[i_tr]>80){
	
	if (Tr_ifStartedAtTarget[i_tr]==1) h_gamma_ene->Fill(genTrE[i_tr]);

      }
    }
  }

  TFile *f_out = new TFile("histos.root","recreate");
  f_out->cd();
  h_ep->Write();
  h_ep_ori->Write();
  h_gamma_ori->Write();
  h_gamma_ene->Write();
}


string IntToString (int a)
{
  ostringstream temp;
  temp<<a;
  return temp.str();
}
