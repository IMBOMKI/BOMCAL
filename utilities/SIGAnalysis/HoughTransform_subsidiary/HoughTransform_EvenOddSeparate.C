#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
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
#include <stdio.h>
#include <time.h>

Double_t ConfTransX(Double_t x, Double_t y);
Double_t ConfTransY(Double_t x, Double_t y);
Double_t HoughTrans(Double_t x, Double_t y, Double_t theta);
Bool_t ifInsideDisk(Double_t x, Double_t y);
std::vector<std::pair<Double_t,Double_t> > makeOrigins(Double_t rad_uncertainty, std::pair<Double_t,Double_t> ref);
Int_t findMaxpoint(std::vector<Int_t> vec);
Bool_t ifInsideBand(Double_t x, Double_t y, Double_t cX, Double_t cY, Double_t outerR, Double_t interR);
Bool_t ifInsideVec(Int_t element, std::vector<Int_t> vec);
Double_t twoPtDistance(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
bool sortFunction (int i,int j) { return (i<j); }
void MakeCluster(std::vector<Int_t> &WireIds, std::vector<std::vector< Int_t> > &ClusterSet, Int_t WireNumberInLayer);
bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size);
std::pair<Double_t, Double_t> getRotatedXY(std::pair<Double_t, Double_t> xy, Double_t deg);
Bool_t ifCircleIsPassing(Double_t rad, Double_t cX, Double_t cY, std::pair<Double_t, Double_t> ref1, std::pair<Double_t, Double_t> ref2);

void HoughTransform_EvenOddSeparate(){
  
  time_t start,end;
  time (&start);

  std::string dir = "../";
  std::string fileName = "trig_ep92_onlyPrimary.root";
  TFile *f = TFile::Open(TString(dir+fileName));
  TTree *t = (TTree*)f->Get("trdata");
  std::string outputdir = "FindedEvents/";
  std::string outputfileName = "find_"+fileName;
  TFile *f_out = new TFile(TString(outputdir+outputfileName),"recreate");
  TTree *t_out = t->CloneTree(0);

  Int_t eventId;
  Int_t nCDCHit;
  Double_t CDCHitX[30000];
  Double_t CDCHitY[30000];
  Double_t CDCHitZ[30000];
  Double_t CDCHitT[30000];
  Double_t CDCEDep[30000];
  Int_t nCALCDCHit;
  Double_t CDCDriftDist[30000];
  Int_t CDCCharge[30000];
  Double_t WireEnd0X[30000];
  Double_t WireEnd0Y[30000];  
  Double_t WireEnd0Z[30000];
  Double_t WireEnd1X[30000];
  Double_t WireEnd1Y[30000];  
  Double_t WireEnd1Z[30000];
  Int_t WireLayerId[30000];
  Int_t WireMaxLayerId;
  Int_t WireId[30000];

  Int_t PDGNumber;
  Double_t genTrE;
  Double_t genTrPx;
  Double_t genTrPy;
  Double_t genTrPz;
  Double_t genTrT;
  Int_t nTrigSTL;
  Int_t nTrigCRK;
  Int_t TrigSTL[1000];
  Int_t TrigCRK[1000];
  Double_t TrigSTLTime[1000];
  Double_t TrigCRKTime[1000];
  Int_t MaxCDCLayerId;
  Int_t nCDCHitPrimary;
  Bool_t ifSingleTurn;
  Bool_t ifMultiTurn;

  Int_t nSTLTrigIndex;
  Int_t STLTrigIndex[200];
  Int_t nCRKTrigIndex;
  Int_t CRKTrigIndex[200];
  //std::vector < std::tuple<std::vector <Int_t>, std::vector <Int_t>, Double_t > > CTHTrigInfo;

  t->SetBranchAddress("eventId", &eventId);
  t->SetBranchAddress("nCDCHit", &nCDCHit);
  t->SetBranchAddress("CDCHitX", CDCHitX);  
  t->SetBranchAddress("CDCHitY", CDCHitY);  
  t->SetBranchAddress("CDCHitZ", CDCHitZ); 
  t->SetBranchAddress("CDCEDep", CDCEDep);   
  t->SetBranchAddress("nCALCDCHit", &nCALCDCHit);
  t->SetBranchAddress("CDCDriftDist", CDCDriftDist);
  t->SetBranchAddress("CDCCharge", CDCCharge);
  t->SetBranchAddress("WireEnd0X", WireEnd0X);
  t->SetBranchAddress("WireEnd0Y", WireEnd0Y);
  t->SetBranchAddress("WireEnd0Z", WireEnd0Z);
  t->SetBranchAddress("WireEnd1X", WireEnd1X);
  t->SetBranchAddress("WireEnd1Y", WireEnd1Y);
  t->SetBranchAddress("WireEnd1Z", WireEnd1Z);
  t->SetBranchAddress("WireLayerId", WireLayerId);
  t->SetBranchAddress("WireMaxLayerId", &WireMaxLayerId);
  t->SetBranchAddress("WireId", WireId);

  t->SetBranchAddress("PDGNumber", &PDGNumber);
  t->SetBranchAddress("genTrE", &genTrE);
  t->SetBranchAddress("genTrPx", &genTrPx);
  t->SetBranchAddress("genTrPy", &genTrPy);
  t->SetBranchAddress("genTrPz", &genTrPz);
  t->SetBranchAddress("genTrT", &genTrT);

  t->SetBranchAddress("nTrigSTL", &nTrigSTL);
  t->SetBranchAddress("nTrigCRK", &nTrigCRK);    
  t->SetBranchAddress("TrigSTL", TrigSTL);
  t->SetBranchAddress("TrigCRK", TrigCRK);
  t->SetBranchAddress("TrigSTLTime", TrigSTLTime);
  t->SetBranchAddress("TrigCRKTime", TrigCRKTime);
  t->SetBranchAddress("MaxCDCLayerId", &MaxCDCLayerId);
  t->SetBranchAddress("nCDCHitPrimary", &nCDCHitPrimary);
  t->SetBranchAddress("ifSingleTurn", &ifSingleTurn);
  t->SetBranchAddress("ifMultiTurn", &ifMultiTurn);

  t->SetBranchAddress("nSTLTrigIndex", &nSTLTrigIndex);
  t->SetBranchAddress("STLTrigIndex", STLTrigIndex);
  t->SetBranchAddress("nCRKTrigIndex", &nCRKTrigIndex);
  t->SetBranchAddress("CRKTrigIndex", CRKTrigIndex);

  //t->SetBranchAddress("CTHTrigInfo", &CTHTrigInfo);
  
  Int_t nRecoHit;
  Double_t RecoWireEnd0X[30000];
  Double_t RecoWireEnd0Y[30000];
  Double_t RecoWireEnd0Z[30000];
  Double_t RecoWireEnd1X[30000];
  Double_t RecoWireEnd1Y[30000];
  Double_t RecoWireEnd1Z[30000];
  Double_t RecoCDCDriftDist[30000];
  Int_t RecoWireLayerId[30000];
  Int_t RecoWireId[30000];
  Int_t RecoCharge;
  Int_t RecoMaxWireLayerId;
  Bool_t Reco_ifCL3;
  Bool_t Reco_ifSingleTurn;
  Bool_t Reco_ifMultiTurn;
  Double_t fittedR_even;
  Double_t fittedR_odd;
  Double_t drEvenToOdd;
  Double_t TruthZ1;

  t_out->Branch("nRecoHit", &nRecoHit, "nRecoHit/I");
  t_out->Branch("RecoWireEnd0X", RecoWireEnd0X, "RecoWireEnd0X[nRecoHit]/D");
  t_out->Branch("RecoWireEnd0Y", RecoWireEnd0Y, "RecoWireEnd0Y[nRecoHit]/D");
  t_out->Branch("RecoWireEnd0Z", RecoWireEnd0Z, "RecoWireEnd0Z[nRecoHit]/D");
  t_out->Branch("RecoWireEnd1X", RecoWireEnd1X, "RecoWireEnd1X[nRecoHit]/D");
  t_out->Branch("RecoWireEnd1Y", RecoWireEnd1Y, "RecoWireEnd1Y[nRecoHit]/D");
  t_out->Branch("RecoWireEnd1Z", RecoWireEnd0Z, "RecoWireEnd1Z[nRecoHit]/D");
  t_out->Branch("RecoCDCDriftDist", RecoCDCDriftDist, "RecoCDCDriftDist[nRecoHit]/D");
  t_out->Branch("RecoWireLayerId", RecoWireLayerId, "RecoWireLayerId[nRecoHit]/I");
  t_out->Branch("RecoWireId", RecoWireId, "RecoWireId[nRecoHit]/I");
  t_out->Branch("RecoCharge", &RecoCharge, "RecoCharge/I");
  t_out->Branch("TruthZ1", &TruthZ1, "TruthZ1/D");
  t_out->Branch("drEvenToOdd", &drEvenToOdd, "drEvenToOdd/D");
  t_out->Branch("RecoMaxWireLayerId", &RecoMaxWireLayerId, "RecoMaxWireLayerId/I");
  t_out->Branch("Reco_ifCL3", &Reco_ifCL3, "Reco_ifCL3/O");
  t_out->Branch("Reco_ifSingleTurn", &Reco_ifSingleTurn, "Reco_ifSingleTurn/O");
  t_out->Branch("Reco_ifMultiTurn", &Reco_ifMultiTurn, "Reco_ifMultiTurn/O");  
  t_out->Branch("fittedR_even", &fittedR_even, "fittedR_even/D");  
  t_out->Branch("fittedR_odd", &fittedR_odd, "fittedR_odd/D");

  TCanvas *c_div = new TCanvas ("c_div", "c divided", 600,900);
  c_div->Divide(2,4);

  //TH1F *maxvote_even_hist = new TH1F("maxvote_even_hist","max_even_vote_hist", 100, 0, 100);
  //TH1F *maxvote_odd_hist = new TH1F("maxvote_odd_hist","max_odd_vote_hist", 100, 0, 100);
  TH1F *fitpT_hist  = new TH1F("fitpT","fitpT",100,60,160) ;
  TH1F *truthpT_hist = new TH1F("truthpT","truthpT",100,60,160) ;
  TH1F *diff_hist = new TH1F("diff","fitpT-truthpT",200,-100,100) ;
  TH1F *rad_even_hist = new TH1F("rad_even_hist","rad_even_hist", 100, 0, 50);
  TH1F *rad_odd_hist = new TH1F("rad_odd_hist","rad_odd_hist", 100, 0, 50);
  TH2F *ref_even_dist = new TH2F("ref_even_dist", "ref_even_dist", 10, -10, 10, 10, -10, 10);
  TH2F *ref_odd_dist = new TH2F("ref_odd_dist", "ref_odd_dist", 10, -10, 10, 10, -10, 10);
 
  TCanvas *c_hits = new TCanvas("c_hits", "c_hits", 1000,1000);
  /*
  TCanvas *c_useful = new TCanvas("c_useful", "c_useful", 1000,1500);
  c_useful->Divide(2,3);
  */
  Int_t NumOfLayers=18;
  Int_t NumOfWiresPerLayer[18]={198,204,210,216,222,228,234,240,246,252,258,264,270,276,282,288,294,300}; // Count only "Actual" Sense Wires
  Double_t LayerRadius[18]={53.0, 54.6, 56.2, 57.8, 59.4, 61.0, 62.6, 64.2, 65.8, 67.4, 69.0, 70.6, 72.2, 73.8, 75.4, 77.0, 78.6, 80.2};
  Double_t LayerStereo_TDR[18]={-67.899, 67.640, -67.384, 67.129, -66.876, 66.625, -66.376, 66.129, -65.884, 65.640, -65.398, 65.158, -64.920, 64.683, -75.132, 74.862, -74.593, 74.326}; // After ICEDUST modified we should use this one!
  Double_t LayerStereo_CAL[18]={-67.004, 66.778, -66.553, 66.327, -66.102, 65.877, -65.652, 65.428, -65.205, 64.982, -64.761, 64.540, -64.319, 64.100, -74.472, 74.220, -73.969, 73.719}; // Up to now, we use this one...
  Double_t ZOffset[18]={-73.6843, -73.9348,-74.2353,-74.5358, -74.7863, -75.0868, 0,0,0,0,0,0,0,0,0,0,0,0};

  //////// Trigger Hodoscope Variables  ////////

  Int_t NumOfCTH=64;
  Double_t STLRad=48.28;
  Double_t STLWidth=9.;
  Double_t STLHeight=0.5;
  Double_t STLTiltAngle=13.;
  Double_t CRKRad=44.78;
  Double_t CRKWidth=9.;
  Double_t CRKHeight=1.;
  Double_t CRKTiltAngle=20.;

  //////// HoughTransform Variables //////// 

  Int_t niter=3;
  Int_t nBins=100;
  Double_t nPt=300.0;   //Double_t nPt=1000.0;
  Double_t rhomax=0.02;
  Double_t rhomin=-0.02;

  //////// Bandwidth Variables ////////

  Int_t TotalHits_Single=0;
  Int_t TotalHits_Multi=0;
  Int_t RecoHits_Single=0;
  Int_t RecoHits_Multi=0;
  
  Double_t bandwidth=5;
  Int_t nIter_band = 10;
  Double_t RecoRate_Single[10];
  Double_t RecoRate_Multi[10];

  /////////////////////////////////////

  for (Int_t i_evt=1; i_evt<7; i_evt++){
    t->GetEntry(i_evt);    

    ///////////////////
    // Pre Track Cut //
    ///////////////////
    
    Bool_t ifThereisOddhit=0;
    Bool_t ifThereisEvenhit=0;
    if (nCALCDCHit<2) continue; // Require At least more than 2 hits
    for (Int_t i_layer=0; i_layer<nCALCDCHit; i_layer++){  
      if ((WireLayerId[i_layer]+1)%2==0) {ifThereisEvenhit=1;}
      else if ((WireLayerId[i_layer]+1)%2==1) {ifThereisOddhit=1;}
    }

    if (ifThereisEvenhit==0 || ifThereisOddhit==0) continue; // Require At least one hit at odd and even layer

    ///////////////////
    // Nullification //
    ///////////////////

    nRecoHit=0;
    memset (RecoWireEnd0X, 0, sizeof(RecoWireEnd0X)); 
    memset (RecoWireEnd0Y, 0, sizeof(RecoWireEnd0Y)); 
    memset (RecoWireEnd0Z, 0, sizeof(RecoWireEnd0Z)); 
    memset (RecoWireEnd1X, 0, sizeof(RecoWireEnd1X)); 
    memset (RecoWireEnd1Y, 0, sizeof(RecoWireEnd1Y)); 
    memset (RecoWireEnd1Z, 0, sizeof(RecoWireEnd1Z)); 
    memset (RecoCDCDriftDist, 0, sizeof(RecoCDCDriftDist));
    memset (RecoWireLayerId, 0, sizeof(RecoWireLayerId));
    RecoCharge=0;
    TruthZ1=0;
    RecoMaxWireLayerId=0;
    Reco_ifCL3=0;
    Reco_ifSingleTurn=0;
    Reco_ifMultiTurn=0;

    
    std::vector<std::pair<Double_t,Double_t> > WireEnd0;
    std::vector<Int_t> LayerId;
    Int_t nEvenhits=0;
    Int_t nOddhits=0;    
    Double_t WireEnd0X_even[30000];
    Double_t WireEnd0Y_even[30000];  
    Double_t WireEnd0X_odd[30000];
    Double_t WireEnd0Y_odd[30000];

    Double_t ConfX[30000];
    Double_t ConfY[30000];
        
    Int_t nConfEven=0;
    Int_t nConfOdd=0;
    Double_t ConfX_even[3000];
    Double_t ConfX_odd[3000];
    Double_t ConfY_even[3000];
    Double_t ConfY_odd[3000];
    
    Int_t deg_index;    
    Int_t rho_index;
    Int_t do_not_use;

    Int_t deg_even_index; 
    Int_t rho_even_index;
    std::pair <Double_t,Double_t> ref_even;
    Int_t deg_odd_index;
    Int_t rho_odd_index;
    std::pair <Double_t,Double_t> ref_odd;

    //////////////////////////////////////////////////////////
    // Sort Even, Odd Hit for Visulaization  &&  Clustering //
    //////////////////////////////////////////////////////////

    std::vector<Int_t> WireIdsPerLayer[18];
    std::vector< std::vector <Int_t> > tmpClusterSet;
    std::vector< std::vector <Int_t> > ClusterSet;

    for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){

      /*
      Double_t theta0=TMath::ATan(WireEnd0Y[i_hit]/WireEnd0X[i_hit]);
      Double_t theta1=TMath::ATan(WireEnd1Y[i_hit]/WireEnd1X[i_hit]);;

      Int_t layernumber=WireLayerId[i_hit];
      Double_t tiltAngle=TMath::ATan(2*LayerRadius[layernumber]/(WireEnd1Z[i_hit]-WireEnd0Z[i_hit]) * TMath::Sin((theta1-theta0)/2));
      //std::cout << layernumber <<"   " << tiltAngle<< "   " << TMath::Tan(LayerStereo[layernumber]/1000) << "   " << theta1-theta0 << std::endl;
      //std::cout << layernumber <<"   " << LayerRadius[layernumber]/(WireEnd1Z[i_hit]-WireEnd0Z[i_hit])* (theta1-theta0) << "   " << TMath::Tan(LayerStereo[layernumber]/1000) << std::endl;

    
      
      Double_t distance = sqrt(pow(WireEnd1X[i_hit]-WireEnd0X[i_hit],2)+pow(WireEnd1Y[i_hit]-WireEnd0Y[i_hit],2)+pow(WireEnd1Z[i_hit]-WireEnd0Z[i_hit],2));
      std::cout <<WireLayerId[i_hit]  <<"   " <<  TMath::ACos((WireEnd1Z[i_hit]-WireEnd0Z[i_hit])/(distance)) << "    " << LayerStereo[WireLayerId[i_hit]]/1000 << "   "<< sqrt(pow(WireEnd0X[i_hit],2)+pow(WireEnd0Y[i_hit],2))  << "   " << LayerRadius[WireLayerId[i_hit]] <<std::endl;
      
      std::cout <<WireLayerId[i_hit]  <<"   " <<  WireEnd0Z[i_hit] << "   " << WireEnd1Z[i_hit] << std::endl;
      */

      // Sort Even, Odd Hit
      if ((WireLayerId[i_hit]+1)%2==0){
	WireEnd0X_even[nEvenhits]=WireEnd0X[i_hit];
	WireEnd0Y_even[nEvenhits]=WireEnd0Y[i_hit];
	nEvenhits++;
      }
      else if((WireLayerId[i_hit]+1)%2==1){
	WireEnd0X_odd[nOddhits]=WireEnd0X[i_hit];
	WireEnd0Y_odd[nOddhits]=WireEnd0Y[i_hit];
	nOddhits++;
      }
      
      // Sort Layer of Hits
      Int_t layerid=WireLayerId[i_hit];
      if(std::find(WireIdsPerLayer[layerid].begin(), WireIdsPerLayer[layerid].end(), WireId[i_hit]) != WireIdsPerLayer[layerid].end()){
	continue;
      }
      WireIdsPerLayer[layerid].push_back(WireId[i_hit]);    	        
    }

    // Sort one more time from smallest to biggest 
    // + Make Clusters and push into ClusterSet
    
    for (Int_t i_layer=0; i_layer<NumOfLayers; i_layer++){
      std::sort (WireIdsPerLayer[i_layer].begin(), WireIdsPerLayer[i_layer].end(), sortFunction);
      for (Int_t i_ele=0; i_ele<WireIdsPerLayer[i_layer].size(); i_ele++){
      }
      MakeCluster(WireIdsPerLayer[i_layer], tmpClusterSet, NumOfWiresPerLayer[i_layer]);
    }

    // Reduce size 1 Cluster
    for (Int_t i_clu=0; i_clu<tmpClusterSet.size(); i_clu++){
      if (tmpClusterSet[i_clu].size()>1){
	ClusterSet.push_back(tmpClusterSet[i_clu]);
      }
    }

   
    ////////////////////////
    //    Main Analysis   //
    ////////////////////////
       
      for (Int_t is_even=0; is_even<2; is_even++){
	  
	Double_t rad_uncertainty=5;
	std::pair <Double_t,Double_t> ref = std::make_pair(0,0);

	for (Int_t it2=0; it2<niter; it2++){
	  
	  std::vector<std::pair<Double_t,Double_t> > origins = makeOrigins(rad_uncertainty,ref);     
	  std::vector<Int_t> vote_max;
	  std::vector<std::pair<Int_t,Int_t> > index_vec; 
	  
	  for (Int_t it=0; it<origins.size(); it++){
	    
	    Double_t oriX=origins[it].first; 
	    Double_t oriY=origins[it].second;
	    
	    if  (ifInsideDisk(oriX,oriY)==0){
	      vote_max.push_back(0);
	    }
	    
            else if (ifInsideDisk(oriX,oriY)==1){
	      
	      /*------------------------------------------
		|                                        |
		|      Conformal Transformation          |
		|                                        |
		-----------------------------------------*/
	      
	      memset(ConfX,0,sizeof(ConfX));
	      memset(ConfY,0,sizeof(ConfY));	    

	      /*
	      memset(ConfX_even,0,sizeof(ConfX_even));
	      memset(ConfY_even,0,sizeof(ConfY_even));
	      memset(ConfX_odd,0,sizeof(ConfX_odd));
	      memset(ConfY_odd,0,sizeof(ConfY_odd));	      
	      nConfEven=0;
	      nConfOdd=0;
	      */
	      for (Int_t i_hit=0; i_hit<nCALCDCHit ; i_hit++){   
		ConfX[i_hit] = ConfTransX(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		ConfY[i_hit] = ConfTransY(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		/*
		if ((WireLayerId[i_hit]+1)%2 == 0){ //even
		  ConfX_even[nConfEven] =  ConfTransX(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		  ConfY_even[nConfEven] =  ConfTransY(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		  nConfEven++;
		}
		else if ((WireLayerId[i_hit]+1)%2 == 1){ //even
		  ConfX_odd[nConfOdd] =  ConfTransX(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		  ConfY_odd[nConfOdd] =  ConfTransY(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		  nConfOdd++;
		}
		*/
	      }

	      
	      /*-------------------------------------------
		|                                         |
		|      Hough Transformation & Voting      |
		|                                         |
		------------------------------------------*/
	      
	      Bool_t if_already_vote[nBins][nBins][2];
	      //TNtuple * hough = new TNtuple("hough","hough", "rho:theta:hitId:wireId:is_even");
	      TH2F *vote = new TH2F("vote_even", "vote_even", nBins, 0, nBins, nBins, 0, nBins);

	      
	      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
		memset(if_already_vote,0,sizeof(if_already_vote));
		for (Int_t i_Pt=0 ; i_Pt<nPt; i_Pt++){
		  
		  Double_t deg = i_Pt/nPt*180.0;
		  Double_t rho =HoughTrans(ConfX[i_hit],ConfY[i_hit],deg);	  
		  Int_t tmp_deg_index = (deg)*nBins/180;
		  Int_t tmp_rho_index = (rho+rhomax)*nBins/(rhomax-rhomin);
	      
		  if (is_even==WireLayerId[i_hit]%2 && if_already_vote[tmp_deg_index][tmp_rho_index][0]==0){
		    if_already_vote[tmp_deg_index][tmp_rho_index][0]=1;
		    vote->Fill(tmp_deg_index,tmp_rho_index);
		  }
		}
	      }
   
	      vote_max.push_back(vote->GetMaximum());	      
	      vote->GetMaximumBin(deg_index, rho_index, do_not_use);
	      index_vec.push_back(std::make_pair(deg_index,rho_index));
	      delete vote;
	    }
	  }

	  Int_t maxIndex=findMaxpoint(vote_max);
	  ref.first = origins[maxIndex].first; 
	  ref.second = origins[maxIndex].second;

	  if (is_even==1){	  
	    ref_even.first = ref.first; 
	    ref_even.second = ref.second;
	    deg_even_index = index_vec[maxIndex].first; 
	    rho_even_index = index_vec[maxIndex].second;

	    if (it2==niter-1){
	      ref_even_dist->Fill(ref_even.first,ref_even.second);}

	  }

	  else if (is_even==0){
	    ref_odd.first = ref.first; 
	    ref_odd.second = ref.second;
	    deg_odd_index = index_vec[maxIndex].first; 
	    rho_odd_index = index_vec[maxIndex].second;

	    if (it2==niter-1){
	      ref_odd_dist->Fill(ref_odd.first,ref_odd.second);}

	  }
	  
	  rad_uncertainty = rad_uncertainty/2.0;
	
	}
      }
      
      Double_t deg_even_eval = 180.0/nBins*(2*deg_even_index+1)/2;
      Double_t rho_even_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_even_index-1)/2;   
      Double_t cX_even = TMath::Cos(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      Double_t cY_even = TMath::Sin(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      
      Double_t deg_odd_eval = 180.0/nBins*(2*deg_odd_index+1)/2;
      Double_t rho_odd_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_odd_index-1)/2;
      Double_t cX_odd = TMath::Cos(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);
      Double_t cY_odd = TMath::Sin(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);

      Double_t rad_even = sqrt(pow(cX_even,2)+pow(cY_even,2));
      Double_t rad_odd = sqrt(pow(cX_odd,2)+pow(cY_odd,2));
      rad_even_hist->Fill(rad_even);
      rad_odd_hist->Fill(rad_odd);
      fittedR_even=rad_even;
      fittedR_odd=rad_odd;


      Double_t fitpT = (rad_even/0.3356+rad_odd/0.3356)/2;
      Double_t truthpT = sqrt(pow(genTrPy,2)+pow(genTrPz,2));
      fitpT_hist->Fill(fitpT);
      truthpT_hist->Fill(truthpT);
      diff_hist->Fill(fitpT-truthpT);

      
      /*-------------------------------------------------
	|                                               |
	|    Count Signal Hit Number in Circle Band     |
	|                                               |
	------------------------------------------------*/
      
      Int_t nRecoHit_even=0;
      Int_t nRecoHit_odd=0;
      Double_t RecoWireEnd0X_even[30000];
      Double_t RecoWireEnd0Y_even[30000];
      Double_t RecoWireEnd0X_odd[30000];
      Double_t RecoWireEnd0Y_odd[30000];
      Double_t outerR_even = rad_even + bandwidth/2.0;
      Double_t interR_even = rad_even - bandwidth/2.0;
      Double_t outerR_odd = rad_odd + bandwidth/2.0;
      Double_t interR_odd = rad_odd - bandwidth/2.0;
      
      RecoMaxWireLayerId=0;
      
      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	if ((WireLayerId[i_hit]+1)%2 == 0){ // even
	  if (ifInsideBand(WireEnd0X[i_hit],WireEnd0Y[i_hit],cX_even+ref_even.first,cY_even+ref_even.second,outerR_even,interR_even)==1){
	    RecoWireEnd0X_even[nRecoHit_even]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_even[nRecoHit_even]=WireEnd0Y[i_hit];
	    nRecoHit_even++;	   

	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd0Z[nRecoHit]=WireEnd0Z[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoWireEnd1Z[nRecoHit]=WireEnd1Z[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];

	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;
	  }
	}
	else if ((WireLayerId[i_hit]+1)%2 == 1){ // odd
	  if (ifInsideBand(WireEnd0X[i_hit],WireEnd0Y[i_hit],cX_odd+ref_odd.first,cY_odd+ref_odd.second,outerR_odd,interR_odd)==1){
	    RecoWireEnd0X_odd[nRecoHit_odd]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_odd[nRecoHit_odd]=WireEnd0Y[i_hit];
	    nRecoHit_odd++;

	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd0Z[nRecoHit]=WireEnd0Z[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoWireEnd1Z[nRecoHit]=WireEnd1Z[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];

	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;
	  }
	}		
      }

      /////////////////////////////////////////////////////
      if (ifSingleTurn==1){
	TotalHits_Single+=nCALCDCHit;
	RecoHits_Single+=(nRecoHit_even+nRecoHit_odd);
      }

      else if (ifMultiTurn==1){
	TotalHits_Multi+=nCALCDCHit;
	RecoHits_Multi+=(nRecoHit_even+nRecoHit_odd);
      }
      ////////////////////////////////////////////////////

      /*-----------------------------------------------------------
	|                                                         |
	|           Count Non-Recognized Hits in Cluster          |
	|                                                         |
	----------------------------------------------------------*/

      std::vector <Int_t> WireIdsToBeReco;

      for (Int_t i_clu=0; i_clu<ClusterSet.size(); i_clu++){
	Bool_t isItRecoCluster=0;

	for(Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){	  
	  if (ifInsideArray(ClusterSet[i_clu][i_ele], RecoWireId, nRecoHit)==1){
	    isItRecoCluster=1;
	    break;
	  }
	}
	
	if (isItRecoCluster==1){
	  for(Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){	  
	    if (ifInsideArray(ClusterSet[i_clu][i_ele], RecoWireId, nRecoHit)==0){
	      WireIdsToBeReco.push_back(ClusterSet[i_clu][i_ele]);
	    }
	  }
	}
      }
  
      //std::cout << "Num Of Ids to be Recognized: " << WireIdsToBeReco.size() << std::endl;

      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	if (ifInsideVec(WireId[i_hit],WireIdsToBeReco)==1){

	  if ((WireLayerId[i_hit]+1)%2 == 0){ // even
	    RecoWireEnd0X_even[nRecoHit_even]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_even[nRecoHit_even]=WireEnd0Y[i_hit];
	    nRecoHit_even++;	   
	    
	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd0Z[nRecoHit]=WireEnd0Z[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoWireEnd1Z[nRecoHit]=WireEnd1Z[i_hit];

	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];
	    
	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;
	    
	  }
	  else if ((WireLayerId[i_hit]+1)%2 == 1){ // odd
	    RecoWireEnd0X_odd[nRecoHit_odd]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_odd[nRecoHit_odd]=WireEnd0Y[i_hit];
	    nRecoHit_odd++;
	    
	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd0Z[nRecoHit]=WireEnd0Z[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoWireEnd1Z[nRecoHit]=WireEnd1Z[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];
	    
	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;	    
	  }			  
	} 
      }
      

      /*--------------------
	|                  |                                      
	|    CL3 Check     |
	|                  |
	-------------------*/

      Double_t abs_cX_even=cX_even+ref_even.first;
      Double_t abs_cY_even=cY_even+ref_even.second;
      Double_t abs_cX_odd=cX_odd+ref_odd.first;
      Double_t abs_cY_odd=cY_odd+ref_odd.second;      
      std::vector< Int_t > domain1_layer;
      std::vector< TVector3 > domain1_wireend0;
      std::vector< TVector3 > domain1_wireend1;
      std::vector< Int_t > domain2_layer;
      std::vector< TVector3 > domain2_wireend0;
      std::vector< TVector3 > domain2_wireend1;


      Double_t slope_even = cY_even/cX_even;
      Double_t slope_odd  = cY_odd/cX_odd;

      for (Int_t i_reco=0; i_reco<nRecoHit; i_reco++){
	TVector3 wireend0(RecoWireEnd0X[i_reco],RecoWireEnd0Y[i_reco],RecoWireEnd0Z[i_reco]);
	TVector3 wireend1(RecoWireEnd1X[i_reco],RecoWireEnd1Y[i_reco],RecoWireEnd1Z[i_reco]);


	if ((RecoWireLayerId[i_reco]+1)%2 == 0){ //even

	  Double_t cross=RecoWireEnd0X[i_reco]*abs_cY_even-RecoWireEnd0Y[i_reco]*abs_cX_even;
	  Double_t LHS= RecoWireEnd0Y[i_reco];
	  Double_t RHS= slope_even*(RecoWireEnd0X[i_reco]-abs_cX_even)+abs_cY_even;
	  if (cross<=0){
	    domain1_layer.push_back(RecoWireLayerId[i_reco]);
	    domain1_wireend0.push_back(wireend0);
	    domain1_wireend1.push_back(wireend1);
	  }
	  else if (cross>0){
	    domain2_layer.push_back(RecoWireLayerId[i_reco]);
	    domain2_wireend0.push_back(wireend0);
	    domain2_wireend1.push_back(wireend1);
	  }
	}

	else if ((RecoWireLayerId[i_reco]+1)%2 == 1){ //odd
	  Double_t cross=RecoWireEnd0X[i_reco]*abs_cY_even-RecoWireEnd0Y[i_reco]*abs_cX_even;
	  Double_t LHS= RecoWireEnd0Y[i_reco];
	  Double_t RHS= slope_odd*(RecoWireEnd0X[i_reco]-abs_cX_odd)+abs_cY_odd;
	  if (cross<=0){
	    domain1_layer.push_back(RecoWireLayerId[i_reco]);
	    domain1_wireend0.push_back(wireend0);
	    domain1_wireend1.push_back(wireend1);
	  }
	  else if (cross>0){
	    domain2_layer.push_back(RecoWireLayerId[i_reco]);
	    domain2_wireend0.push_back(wireend0);
	    domain2_wireend1.push_back(wireend1);
	  }
	}
      }

      if (ifInsideVec(0,domain1_layer)==1 && ifInsideVec(1,domain1_layer)==1 && ifInsideVec(2,domain1_layer)==1){
	if (ifInsideVec(0,domain2_layer)==1 && ifInsideVec(1,domain2_layer)==1 && ifInsideVec(2,domain2_layer)==1){
	  Reco_ifCL3=1;
	}
      }

      /*--------------------------------
	|                              |
	|    Charge Identification     |
	|                              |
	-------------------------------*/
     
      std::vector < Int_t > STLUpIndex;
      std::vector < Int_t > STLDownIndex;
      std::vector < Int_t > CRKUpIndex;
      std::vector < Int_t > CRKDownIndex;  // Up & Down Indices Normalized into NumOfCTH Scale (<64)

      std::vector < Int_t > STLRecoIndex;
      std::vector < Int_t > CRKRecoIndex;  // Indices where Hough Circles are passing through
      std::vector < Bool_t > isIndexClockwise;

      for (Int_t i_stl=0; i_stl<nSTLTrigIndex; i_stl++){      
	Int_t index=STLTrigIndex[i_stl];

	if (index<NumOfCTH){  //Up 
	  STLUpIndex.push_back(index);
	}
	else if (index>=NumOfCTH){ //Down 
	  index=(Int_t(1.5*NumOfCTH)-index%NumOfCTH)%NumOfCTH;
	  STLDownIndex.push_back(index);
	}

	std::pair<Double_t, Double_t> STLx1y1 = std::make_pair(-STLWidth/2.,-STLHeight/2.);
	std::pair<Double_t, Double_t> STLx1y2 = std::make_pair(-STLWidth/2.,STLHeight/2.);
	std::pair<Double_t, Double_t> STLx2y1 = std::make_pair(STLWidth/2.,-STLHeight/2.);
	std::pair<Double_t, Double_t> STLx2y2 = std::make_pair(STLWidth/2.,STLHeight/2.);
	Double_t STLRotAngle=-(90-(index+1/2.)*360./NumOfCTH+STLTiltAngle);
	Double_t xTrans=STLRad*TMath::Cos((index+1/2.)*360./NumOfCTH*TMath::Pi()/180.);
	Double_t yTrans=STLRad*TMath::Sin((index+1/2.)*360./NumOfCTH*TMath::Pi()/180.);

	STLx1y1=getRotatedXY(STLx1y1, STLRotAngle);
	STLx1y2=getRotatedXY(STLx1y2, STLRotAngle);
	STLx2y1=getRotatedXY(STLx2y1, STLRotAngle);
	STLx2y2=getRotatedXY(STLx2y2, STLRotAngle);

	STLx1y1=std::make_pair(STLx1y1.first+xTrans, STLx1y1.second+yTrans);
	STLx1y2=std::make_pair(STLx1y2.first+xTrans, STLx1y2.second+yTrans);
	STLx2y1=std::make_pair(STLx2y1.first+xTrans, STLx2y1.second+yTrans);
	STLx2y2=std::make_pair(STLx2y2.first+xTrans, STLx2y2.second+yTrans);
	
	std::pair <Double_t, Double_t> ref1 = std::make_pair((STLx1y1.first+STLx1y2.first)/2,(STLx1y1.second+STLx1y2.second)/2);
	std::pair <Double_t, Double_t> ref2 = std::make_pair((STLx2y1.first+STLx2y2.first)/2,(STLx2y1.second+STLx2y2.second)/2);
	std::pair <Double_t, Double_t> refMid = std::make_pair((ref1.first+ref2.first)/2,(ref1.second+ref2.second)/2);
	Double_t avg_cX=(abs_cX_even+abs_cX_odd)/2;
	Double_t avg_cY=(abs_cY_even+abs_cY_odd)/2;

	if (ifCircleIsPassing(rad_even, abs_cX_even, abs_cY_even, ref1, ref2)==1 || ifCircleIsPassing(rad_odd, abs_cX_odd, abs_cY_odd, ref1, ref2)==1){
	  STLRecoIndex.push_back(index);
	  Double_t cross = refMid.first*avg_cY-avg_cX*refMid.second;
	  if (cross>=0){
	    isIndexClockwise.push_back(1);
	  }
	  else if (cross<0){
	    isIndexClockwise.push_back(0);
	  }
	}
      }

      for (Int_t i_crk=0; i_crk<nCRKTrigIndex; i_crk++){      
	Int_t index=CRKTrigIndex[i_crk];

	if (index<NumOfCTH){  //Up 
	  CRKUpIndex.push_back(index);
	}
	else if (index>=NumOfCTH){ //Down 
	  index=(Int_t(1.5*NumOfCTH)-index%NumOfCTH)%NumOfCTH;
	  CRKDownIndex.push_back(index);
	}

	std::pair<Double_t, Double_t> CRKx1y1 = std::make_pair(-CRKWidth/2.,-CRKHeight/2.);
	std::pair<Double_t, Double_t> CRKx1y2 = std::make_pair(-CRKWidth/2.,CRKHeight/2.);
	std::pair<Double_t, Double_t> CRKx2y1 = std::make_pair(CRKWidth/2.,-CRKHeight/2.);
	std::pair<Double_t, Double_t> CRKx2y2 = std::make_pair(CRKWidth/2.,CRKHeight/2.);
	Double_t CRKRotAngle=-(90-(index+1/2.)*360./NumOfCTH+CRKTiltAngle);
	Double_t xTrans=CRKRad*TMath::Cos((index+1/2.)*360./NumOfCTH*TMath::Pi()/180.);
	Double_t yTrans=CRKRad*TMath::Sin((index+1/2.)*360./NumOfCTH*TMath::Pi()/180.);

	CRKx1y1=getRotatedXY(CRKx1y1, CRKRotAngle);
	CRKx1y2=getRotatedXY(CRKx1y2, CRKRotAngle);
	CRKx2y1=getRotatedXY(CRKx2y1, CRKRotAngle);
	CRKx2y2=getRotatedXY(CRKx2y2, CRKRotAngle);

	CRKx1y1=std::make_pair(CRKx1y1.first+xTrans, CRKx1y1.second+yTrans);
	CRKx1y2=std::make_pair(CRKx1y2.first+xTrans, CRKx1y2.second+yTrans);
	CRKx2y1=std::make_pair(CRKx2y1.first+xTrans, CRKx2y1.second+yTrans);
	CRKx2y2=std::make_pair(CRKx2y2.first+xTrans, CRKx2y2.second+yTrans);
	
	std::pair <Double_t, Double_t> ref1 = std::make_pair((CRKx1y1.first+CRKx1y2.first)/2,(CRKx1y1.second+CRKx1y2.second)/2);
	std::pair <Double_t, Double_t> ref2 = std::make_pair((CRKx2y1.first+CRKx2y2.first)/2,(CRKx2y1.second+CRKx2y2.second)/2);
	std::pair <Double_t, Double_t> refMid = std::make_pair((ref1.first+ref2.first)/2,(ref1.second+ref2.second)/2);
	Double_t avg_cX=(abs_cX_even+abs_cX_odd)/2;
	Double_t avg_cY=(abs_cY_even+abs_cY_odd)/2;

	if (ifCircleIsPassing(rad_even, abs_cX_even, abs_cY_even, ref1, ref2)==1 || ifCircleIsPassing(rad_odd, abs_cX_odd, abs_cY_odd, ref1, ref2)==1){
	  CRKRecoIndex.push_back(index);
	  Double_t cross = refMid.first*avg_cY-avg_cX*refMid.second;
	  if (cross>=0){
	    isIndexClockwise.push_back(1);
	  }
	  else if (cross<0){
	    isIndexClockwise.push_back(0);
	  }
	}
      }      

      if (isIndexClockwise.empty()==1){
	RecoCharge=0; // Non-Classified
      }
      else{
	if (std::find(isIndexClockwise.begin(), isIndexClockwise.end(), 0) == isIndexClockwise.end()){ //vector only contains 1
	  RecoCharge=1;
	}
	else if (std::find(isIndexClockwise.begin(), isIndexClockwise.end(), 1) == isIndexClockwise.end()){ //vector only contains 0
	  RecoCharge=-1;
	}
	else {
	  RecoCharge=0;
	}
      }      

      /*------------------------------------
	|                                  |
	|    Other Auxiliary Variables     | 
	|                                  |
	-----------------------------------*/

      drEvenToOdd = sqrt(pow(abs_cX_even-abs_cX_odd,2)+pow(abs_cY_even-abs_cY_odd,2));
      TruthZ1=CDCHitZ[0];

      /*-------------------
	|                 |
	|    Printing     | 
	|                 |
	------------------*/

      std::cout << "-------------------------" << std::endl;
      std::cout << i_evt <<"-th Event (EventId: "<< eventId <<")" << std::endl;           
      std::cout << "MC Hit NDF: " << nCALCDCHit << std::endl;
      std::cout << "Reco Charge: " << RecoCharge << std::endl;
      std::cout << "STL Trigger Index: ";
      for (Int_t i_stl=0; i_stl<nSTLTrigIndex ; i_stl++){
	std::cout << STLTrigIndex[i_stl] << " ";
      }
      std::cout << std::endl;
      std::cout << "CRK Trigger Index: ";
      for (Int_t i_crk=0; i_crk<nCRKTrigIndex ; i_crk++){
	std::cout << CRKTrigIndex[i_crk] << " ";
      }
      std::cout << std::endl;

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      //                                                                                               //
      //                                         PLOTTING SECTION                                      //  
      //                                                                                               //
      ///////////////////////////////////////////////////////////////////////////////////////////////////

      /*-----------------------------------------
	|                                         |
	|   Wire Hit Drawing with Hough Circle    |
	|                                         |
	------------------------------------------*/
      
      c_hits->cd();
      c_hits->DrawFrame(-90,-90,90,90);
      
      // Wire Layers
      
      Double_t WireRad[18] = {53.00,54.60,56.20,57.80,59.40,61.00,62.60,64.20,65.80,67.40,69.00,70.60,72.20,73.80,75.40,77.00,78.60,80.20};
      
      for (Int_t i=0; i<18; i++){
	TEllipse *WireCircle = new TEllipse(0,0,WireRad[i],WireRad[i]);
	WireCircle->SetFillColor(0);
	WireCircle->SetFillStyle(4000);
	WireCircle->SetLineColor(40);
	WireCircle->Draw();
      }
      
      // Stopping Targets
      
      Double_t DiskRad=10.00;
      TEllipse *Disk = new TEllipse(0,0,DiskRad,DiskRad);
      Disk->SetFillColor(40);
      Disk->SetLineColor(40);
      Disk->Draw();  

      // Trigger Hodoscopes
            
      for (Int_t i_cth=0; i_cth<NumOfCTH; i_cth++){
	std::pair<Double_t, Double_t> STLx1y1 = std::make_pair(-STLWidth/2.,-STLHeight/2.);
	std::pair<Double_t, Double_t> STLx1y2 = std::make_pair(-STLWidth/2.,STLHeight/2.);
	std::pair<Double_t, Double_t> STLx2y1 = std::make_pair(STLWidth/2.,-STLHeight/2.);
	std::pair<Double_t, Double_t> STLx2y2 = std::make_pair(STLWidth/2.,STLHeight/2.);
	Double_t STLRotAngle=-(90-(i_cth+1/2.)*360./NumOfCTH+STLTiltAngle);
	Double_t xTrans=STLRad*TMath::Cos((i_cth+1/2.)*360./NumOfCTH*TMath::Pi()/180.);
	Double_t yTrans=STLRad*TMath::Sin((i_cth+1/2.)*360./NumOfCTH*TMath::Pi()/180.);

	STLx1y1=getRotatedXY(STLx1y1, STLRotAngle);
	STLx1y2=getRotatedXY(STLx1y2, STLRotAngle);
	STLx2y1=getRotatedXY(STLx2y1, STLRotAngle);
	STLx2y2=getRotatedXY(STLx2y2, STLRotAngle);

	STLx1y1=std::make_pair(STLx1y1.first+xTrans, STLx1y1.second+yTrans);
	STLx1y2=std::make_pair(STLx1y2.first+xTrans, STLx1y2.second+yTrans);
	STLx2y1=std::make_pair(STLx2y1.first+xTrans, STLx2y1.second+yTrans);
	STLx2y2=std::make_pair(STLx2y2.first+xTrans, STLx2y2.second+yTrans);

	Double_t x[5]={STLx1y1.first,  STLx1y2.first,  STLx2y2.first,  STLx2y1.first,  STLx1y1.first};
	Double_t y[5]={STLx1y1.second, STLx1y2.second, STLx2y2.second, STLx2y1.second, STLx1y1.second};

	TPolyLine *STLBox=new TPolyLine(5,x,y);
	STLBox->SetFillColor(0);
	STLBox->SetFillStyle(4000);
		
	if (ifInsideVec(i_cth,STLUpIndex)==1 || ifInsideVec(i_cth,STLDownIndex)==1){
	  STLBox->SetFillColor(8);
	  STLBox->SetFillStyle(1001);
	}
	
	STLBox->Draw("f");
	STLBox->Draw();
      }
      for (Int_t i_cth=0; i_cth<NumOfCTH; i_cth++){
	std::pair<Double_t, Double_t> CRKx1y1 = std::make_pair(-CRKWidth/2.,-CRKHeight/2.);
	std::pair<Double_t, Double_t> CRKx1y2 = std::make_pair(-CRKWidth/2.,CRKHeight/2.);
	std::pair<Double_t, Double_t> CRKx2y1 = std::make_pair(CRKWidth/2.,-CRKHeight/2.);
	std::pair<Double_t, Double_t> CRKx2y2 = std::make_pair(CRKWidth/2.,CRKHeight/2.);
	Double_t CRKRotAngle=-(90-(i_cth+1/2.)*360./NumOfCTH+CRKTiltAngle);
	Double_t xTrans=CRKRad*TMath::Cos((i_cth+1/2.)*360./NumOfCTH*TMath::Pi()/180.);
	Double_t yTrans=CRKRad*TMath::Sin((i_cth+1/2.)*360./NumOfCTH*TMath::Pi()/180.);

	CRKx1y1=getRotatedXY(CRKx1y1, CRKRotAngle);
	CRKx1y2=getRotatedXY(CRKx1y2, CRKRotAngle);
	CRKx2y1=getRotatedXY(CRKx2y1, CRKRotAngle);
	CRKx2y2=getRotatedXY(CRKx2y2, CRKRotAngle);

	CRKx1y1=std::make_pair(CRKx1y1.first+xTrans, CRKx1y1.second+yTrans);
	CRKx1y2=std::make_pair(CRKx1y2.first+xTrans, CRKx1y2.second+yTrans);
	CRKx2y1=std::make_pair(CRKx2y1.first+xTrans, CRKx2y1.second+yTrans);
	CRKx2y2=std::make_pair(CRKx2y2.first+xTrans, CRKx2y2.second+yTrans);

	Double_t x[5]={CRKx1y1.first,  CRKx1y2.first,  CRKx2y2.first,  CRKx2y1.first,  CRKx1y1.first};
	Double_t y[5]={CRKx1y1.second, CRKx1y2.second, CRKx2y2.second, CRKx2y1.second, CRKx1y1.second};

	TPolyLine *CRKBox=new TPolyLine(5,x,y);

	CRKBox->SetFillColor(0);
	CRKBox->SetFillStyle(4000);

	if (ifInsideVec(i_cth,CRKUpIndex)==1 || ifInsideVec(i_cth,CRKDownIndex)==1){
	  CRKBox->SetFillColor(9);
	  CRKBox->SetFillStyle(1001);
	}

	CRKBox->Draw("f");
	CRKBox->Draw();

      }


      // Hough Circles (Even, Odd)
      /*
      TEllipse *circle_even = new TEllipse(abs_cX_even,abs_cY_even,rad_even,rad_even);
      TEllipse *circle_odd = new TEllipse(abs_cX_odd,abs_cY_odd,rad_odd,rad_odd);
      TEllipse *center_even = new TEllipse(abs_cX_even,abs_cY_even,0.1,0.1);
      TEllipse *center_odd = new TEllipse(abs_cX_odd,abs_cY_odd,0.1,0.1);

      circle_even->SetFillColor(0);
      circle_even->SetFillStyle(4000);
      circle_even->SetLineColor(4);
      circle_even->SetLineWidth(1);
      circle_even->Draw();
      
      circle_odd->SetFillColor(0);
      circle_odd->SetFillStyle(4000);
      circle_odd->SetLineColor(2);      
      circle_odd->SetLineWidth(1);
      circle_odd->Draw();
      */
      // Hits
      
      TGraph *grEvenhits = new TGraph(nEvenhits, WireEnd0X_even, WireEnd0Y_even);
      grEvenhits->SetTitle("Evenhits");
      grEvenhits->SetMarkerStyle(20);
      grEvenhits->SetMarkerSize(1);
      grEvenhits->SetMarkerColor(4);  //even - blue
      grEvenhits->Draw("P");
      
      TGraph *grOddhits = new TGraph(nOddhits, WireEnd0X_odd, WireEnd0Y_odd);
      grOddhits->SetTitle("Oddhits");
      grOddhits->SetMarkerStyle(20);
      grOddhits->SetMarkerSize(1);
      grOddhits->SetMarkerColor(2);  //odd - red
      grOddhits->Draw("P");
      
      // Recognized Hits
      
      TGraph *grRecoEvenhits = new TGraph(nRecoHit_even, RecoWireEnd0X_even, RecoWireEnd0Y_even);
      grRecoEvenhits->SetTitle("RecoEvenhits");  
      grRecoEvenhits->SetMarkerStyle(20);
      grRecoEvenhits->SetMarkerSize(1);
      grRecoEvenhits->SetMarkerColor(38); // Reco_even - right blue
      grRecoEvenhits->Draw("P");
      
      TGraph *grRecoOddhits = new TGraph(nRecoHit_odd, RecoWireEnd0X_odd, RecoWireEnd0Y_odd);
      grRecoOddhits->SetTitle("RecoOddhits");
      grRecoOddhits->SetMarkerStyle(20);
      grRecoOddhits->SetMarkerSize(1);
      grRecoOddhits->SetMarkerColor(46); // Reco_odd - right red
      grRecoOddhits->Draw("P");      
      
      /*
      TString figName="hits"+std::to_string(i_evt)+".pdf";
      c_hits->Print(figName);
      */
      
      // Center crossing line
      
      /*
      Int_t nPoints = 8;
      Double_t xLine_even[nPoints];
      Double_t yLine_even[nPoints];
      Double_t xLine_odd[nPoints];
      Double_t yLine_odd[nPoints];

      for (Int_t i_line=0; i_line<nPoints; i_line++) {	
	xLine_even[i_line] = i_line*10;
	yLine_even[i_line] = slope_even*(xLine_even[i_line]-abs_cX_even)+abs_cY_even;
      }
      for (Int_t i_line=0; i_line<nPoints; i_line++) {
	xLine_odd[i_line] = i_line*10;
	yLine_odd[i_line] = slope_odd*(xLine_odd[i_line]-abs_cX_odd)+abs_cY_odd;
      }
      
      
      TGraph *crossingLine_even = new TGraph(nPoints, xLine_even, yLine_even);
      crossingLine_even->SetLineColor(4);
      crossingLine_even->Draw("L");

      TGraph *crossingLine_odd = new TGraph(nPoints, xLine_odd, yLine_odd);
      crossingLine_odd->SetLineColor(2);
      crossingLine_odd->Draw("L");
      */


       /*-----------------------------------
	|                                 |
	|   Conformal Transformed Hits    |
	|                                 |
	----------------------------------*/
      /*
      c_useful->cd(5);

      TGraph *conf_odd = new TGraph(nConfOdd, ConfX_odd, ConfY_odd);
      conf_odd->SetMarkerStyle(20);
      conf_odd->SetMarkerColor(2);
      conf_odd->Draw("AP");

      TGraph *conf_even = new TGraph(nConfEven, ConfX_even, ConfY_even);
      conf_even->GetYaxis()->SetRangeUser(-0.02,0.02);
      conf_even->GetXaxis()->SetRangeUser(-0.02,0.02);
      conf_even->SetMarkerStyle(20);
      conf_even->SetMarkerColor(4);
      conf_even->Draw("P");
      */

      /*----------------------------------------
	|                                      |
	|   Hough Transform  && Voting Plot    |
	|                                      |
	---------------------------------------*/
      /*
      TH2F *vote_plot_even = new TH2F("vote_plot_even", "vote_plot_even", nBins, 0, nBins, nBins, 0, nBins);
      TH2F *vote_plot_odd = new TH2F("vote_plot_odd", "vote_plot_odd", nBins, 0, nBins, nBins, 0, nBins);
      TNtuple * hough = new TNtuple("hough","hough", "rho:theta");
      TNtuple * hough_odd = new TNtuple("hough_odd","hough_odd", "rho:theta");

      Bool_t if_already_vote[nBins][nBins][2];
            
      is_even=1;
      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	memset(if_already_vote,0,sizeof(if_already_vote));
	for (Int_t i_Pt=0 ; i_Pt<nPt; i_Pt++){
	  
	  Double_t deg = i_Pt/nPt*180.0;
	  Double_t rho =HoughTrans(ConfX[i_hit],ConfY[i_hit],deg);	  
	  Int_t tmp_deg_index = (deg)*nBins/180;
	  Int_t tmp_rho_index = (rho+rhomax)*nBins/(rhomax-rhomin);

	  if (is_even==(WireLayerId[i_hit]+1)%2){
	    hough->Fill(rho,deg);
	  }

	  else if (is_even!=(WireLayerId[i_hit]+1)%2){
	    hough_odd->Fill(rho,deg);	  
	  }	  


	  if (is_even==WireLayerId[i_hit]%2 && if_already_vote[tmp_deg_index][tmp_rho_index][0]==0){
	  
	    if_already_vote[tmp_deg_index][tmp_rho_index][0]=1;
	    vote_plot_even->Fill(tmp_deg_index,tmp_rho_index);
	  }
	}
      }
      is_even=0;
      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	memset(if_already_vote,0,sizeof(if_already_vote));
	for (Int_t i_Pt=0 ; i_Pt<nPt; i_Pt++){
	  
	  Double_t deg = i_Pt/nPt*180.0;
	  Double_t rho =HoughTrans(ConfX[i_hit],ConfY[i_hit],deg);	  
	  Int_t tmp_deg_index = (deg)*nBins/180;
	  Int_t tmp_rho_index = (rho+rhomax)*nBins/(rhomax-rhomin);

	  if (is_even==(WireLayerId[i_hit])%2 && if_already_vote[tmp_deg_index][tmp_rho_index][0]==0){
	    if_already_vote[tmp_deg_index][tmp_rho_index][0]=1;
	    vote_plot_odd->Fill(tmp_deg_index,tmp_rho_index);
	  }
	}
      }

      std::cout << "       sadadasdsds  " <<  hough->GetEvent() << std::endl;

      c_useful->cd(1);
      hough->Draw("rho:theta");
      TGraph *hough_gr = new TGraph(hough->GetSelectedRows(), hough->GetV2(), hough->GetV1());
      hough_gr->GetXaxis()->SetTitle("");
      hough_gr->GetXaxis()->SetRangeUser(0,180);
      hough_gr->GetYaxis()->SetTitle("");
      hough_gr->GetYaxis()->SetRangeUser(-0.025,0.025);
      hough_gr->SetLineWidth(1);
      hough_gr->Draw("AP");
      hough->Print();

      c_useful->cd(2);
      hough_odd->Draw("rho:theta");
      TGraph *hough_gr_odd = new TGraph(hough_odd->GetSelectedRows(), hough_odd->GetV2(), hough_odd->GetV1());
      hough_gr_odd->GetXaxis()->SetTitle("");
      hough_gr_odd->GetXaxis()->SetRangeUser(0,180);
      hough_gr_odd->GetYaxis()->SetTitle("");
      hough_gr_odd->GetYaxis()->SetRangeUser(-0.025,0.025);
      hough_gr_odd->SetLineWidth(1);
      hough_gr_odd->Draw("AP");
      hough_odd->Print();

      c_useful->cd(3);    
      vote_plot_even->Draw("colz");
      c_useful->cd(4);
      vote_plot_odd->Draw("colz");
      c_useful->Print("figures.pdf");
      */
      ///////////////////////////////////////////////////////////////////////////////////////////


      t_out->Fill();
  }


  
  c_div->cd(1);
  diff_hist->SetLineColor(4);
  diff_hist->SetFillColor(4);
  diff_hist->SetFillStyle(3004);
  diff_hist->Draw();

  c_div->cd(2);
  rad_even_hist->SetLineColor(4);
  rad_even_hist->SetFillColor(4);
  rad_even_hist->SetFillStyle(3004);
  rad_even_hist->Draw();
  rad_odd_hist->SetLineColor(2);
  rad_odd_hist->SetFillColor(2);
  rad_odd_hist->SetFillStyle(3005);
  rad_odd_hist->Draw("same");

  c_div->cd(3);
  truthpT_hist->SetLineColor(1);
  truthpT_hist->SetFillColor(4);
  truthpT_hist->SetFillStyle(3004);
  truthpT_hist->Draw();

  c_div->cd(4);
  fitpT_hist->SetLineColor(1);
  fitpT_hist->SetFillColor(4);
  fitpT_hist->SetFillStyle(3004);
  fitpT_hist->Draw();

  c_div->cd(5);
  ref_even_dist->Draw("colz");
  c_div->cd(6);
  ref_odd_dist->Draw("colz");
  
  std::cout << "Finish!" << std::endl;
  time (&end);
  double elapsedTime = difftime (end,start);
  std::cout << "ElapsedTime: " << elapsedTime << "sec" << std::endl;

  
  f_out->cd();
  f_out->Write();
  f_out->Close();
  f->cd();
  f->Close();

}

Double_t ConfTransX(Double_t x, Double_t y){
  return x/(pow(x,2)+pow(y,2));
} 

Double_t ConfTransY(Double_t x, Double_t y){
  return y/(pow(x,2)+pow(y,2));
} 

Double_t HoughTrans(Double_t x, Double_t y, Double_t theta){
  return x*TMath::Cos(theta*TMath::Pi()/180)+y*TMath::Sin(theta*TMath::Pi()/180);
}

Bool_t ifInsideDisk(Double_t x, Double_t y){
  if (pow(x,2)+pow(y,2)<100){ return 1;}
  else {return 0;}
};

std::vector<std::pair<Double_t,Double_t> > makeOrigins(Double_t rad_uncertainty, std::pair<Double_t,Double_t> ref){
  std::vector<std::pair<Double_t,Double_t> > origins;
  origins.push_back(std::make_pair(ref.first, ref.second));  
  origins.push_back(std::make_pair(ref.first+rad_uncertainty, ref.second+rad_uncertainty));
  origins.push_back(std::make_pair(ref.first-rad_uncertainty, ref.second+rad_uncertainty));
  origins.push_back(std::make_pair(ref.first-rad_uncertainty, ref.second-rad_uncertainty));
  origins.push_back(std::make_pair(ref.first+rad_uncertainty, ref.second-rad_uncertainty));
  origins.push_back(std::make_pair(ref.first, ref.second+rad_uncertainty));
  origins.push_back(std::make_pair(ref.first, ref.second-rad_uncertainty));
  origins.push_back(std::make_pair(ref.first+rad_uncertainty, ref.second));
  origins.push_back(std::make_pair(ref.first-rad_uncertainty, ref.second));
  return origins;
}

Int_t findMaxpoint(std::vector<Int_t> vec){
  Int_t maxpoint;
  Int_t val=0;
  for (int i=0; i<vec.size(); i++){
    if (vec[i]>val){
      val=vec[i];
      maxpoint=i;
    }
  }
  return maxpoint;
}

Bool_t ifInsideBand(Double_t x, Double_t y, Double_t cX, Double_t cY, Double_t outerR, Double_t interR){
  if ((pow(x-cX,2)+pow(y-cY,2))<pow(outerR,2) && (pow(x-cX,2)+pow(y-cY,2))>pow(interR,2)){
    return 1;
  }
  return 0;
}

Bool_t ifInsideVec(Int_t element, std::vector<Int_t> vec){
  for (Int_t i=0; i<vec.size(); i++){
    if (vec[i]==element){
      return 1;
    }
  }
  return 0;
}

Double_t twoPtDistance(Double_t x1, Double_t y1, Double_t x2, Double_t y2){
  return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
}


void MakeCluster(std::vector<Int_t> &WireIds, std::vector<std::vector< Int_t> > &ClusterSet, Int_t WireNumberInLayer)
{
  if (!WireIds.empty()){
    std::vector<Int_t> tmpCluster;
    std::vector<Int_t> firstCluster;
    tmpCluster.push_back(WireIds[0]);
    firstCluster.push_back(WireIds[0]);

    for(Int_t i_id = 1; i_id < WireIds.size(); ++i_id){
      if(WireIds[i_id] == WireIds[i_id-1] + 1){
	tmpCluster.push_back(WireIds[i_id]);
	
	if (i_id==WireIds.size()-1){ // Last element 
	  // 1. Check if it satisfies periodic boundary condition with first id of layer
	  if (WireIds[i_id]-WireIds[0]+1==WireNumberInLayer){
	    firstCluster.insert(firstCluster.end(), tmpCluster.begin(),tmpCluster.end());	    
	    ClusterSet.push_back(firstCluster);	
	    //	    std::cout << "Found Periodic Boundary!" << std::endl;
	  }
	  // 2. Otherwise, just push it.
	  else{
	    ClusterSet.push_back(firstCluster);
	    ClusterSet.push_back(tmpCluster);
	  }
  
	}

      }
      else if (WireIds[i_id] != WireIds[i_id-1] + 1){
	if (std::find(tmpCluster.begin(), tmpCluster.end(), WireIds[0]) != tmpCluster.end()){
	  firstCluster=tmpCluster;
	}
	else {
	  ClusterSet.push_back(tmpCluster);
	}
	tmpCluster.clear();
	tmpCluster.push_back(WireIds[i_id]);
      }
    }

  }
}

bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size){
  for (Int_t i=0; i<arr_size; i++){
    if (reco_ids[i]==wireid) return 1;
  }
  return 0;
}

std::pair<Double_t, Double_t> getRotatedXY(std::pair<Double_t, Double_t> xy, Double_t deg){
  Double_t RotX = xy.first*TMath::Cos(deg*TMath::Pi()/180)-xy.second*TMath::Sin(deg*TMath::Pi()/180);
  Double_t RotY = xy.second*TMath::Cos(deg*TMath::Pi()/180)+xy.first*TMath::Sin(deg*TMath::Pi()/180);  
  std::pair<Double_t, Double_t> RotatedXY= std::make_pair(RotX, RotY);
  return RotatedXY;
};

Bool_t ifCircleIsPassing(Double_t rad, Double_t cX, Double_t cY, std::pair<Double_t, Double_t> ref1, std::pair<Double_t, Double_t> ref2){
  Double_t dist1 = pow(ref1.first-cX,2)+pow(ref1.second-cY,2);
  Double_t dist2 = pow(ref2.first-cX,2)+pow(ref2.second-cY,2);
  if (dist1<pow(rad,2) && dist2>pow(rad,2)) return 1;
  else if (dist1>pow(rad,2) && dist2<pow(rad,2)) return 1;
  return 0;
};

