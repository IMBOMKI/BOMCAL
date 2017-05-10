#include <IHoughTransform.hxx>

#include "HEPUnits.hxx"
#include <vector>
#include <iostream>
#include <algorithm>

#include <IOADatabase.hxx>
#include <IReconBase.hxx>
#include <IReconTrackCand.hxx>
#include <ITrackState.hxx>

#include <ICOMETEvent.hxx>
#include <ICOMETLog.hxx>

#include <IHitSelection.hxx>
#include <IMCHit.hxx>
#include <IG4HitSegment.hxx>
#include <IG4HitCalo.hxx>
#include <IHit.hxx>
#include <IGeomInfo.hxx>
#include <IGeometryId.hxx>
#include <COMETGeomId.hxx>
#include <IGeomIdManager.hxx>
#include <ICTHGeomId.hxx>
#include <ICTHChannelId.hxx>
#include <IChannelId.hxx>
#include <IGeometryDatabase.hxx>
#include <ICDCGeom.hxx>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TRotation.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2F.h"
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
#include "TEllipse.h"
#include "TPolyLine.h"

bool IntComparison (int i,int j);

void MakeCluster(std::vector<int> &WireIds, std::vector<std::vector< int> > &ClusterSet, int WireNumberInLayer);

bool ifInsideDisk(double x, double y);
bool ifInsideBand(double x, double y, double cX, double cY, double outerR, double interR);
bool ifInsideVec(int element, std::vector<int> vec);
bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size);
Int_t findMaxPoint(std::vector<Int_t> vec);
std::pair<Double_t, Double_t> getRotatedXY(std::pair<Double_t, Double_t> xy, Double_t deg);
Bool_t ifCircleIsPassing(Double_t rad, Double_t cX, Double_t cY, std::pair<Double_t, Double_t> ref1, std::pair<Double_t, Double_t> ref2);

IHoughTransform::IHoughTransform(const char* name = "IHoughTransform", const char* title="houghtransform")
  :ITracking(name,title),
   fnIter(3),
   fnBins(100),
   fnPt(300.0),
   fRhoMax(0.02),
   fRhoMin(-0.02),
   fBandWidth(8),
   fRadUncertainty(5),
   fRef(std::make_pair(0,0)),
   fnRecoHit(0),
   fnRecoHit_even(0),
   fnRecoHit_odd(0),
   fRecoCL3(0),
   fReco2DCharge(0),

   /*** Geometry Variables ***/
   fDiskRad(10.0),
   fCTHSegNum(64),
   fNumOfLayers(18),
   fScintRad(48.28),
   fScintWidth(9.0),
   fScintHeight(0.5),
   fScintTiltAngle(13.),
   fCherenRad(44.78),
   fCherenWidth(9.0),
  fCherenHeight(1.0),
  fCherenTiltAngle(20.)
{ 
  fNumOfWiresPerLayer = new int[18]{198,204,210,216,222,228,234,240,246,252,258,264,270,276,282,288,294,300};
  fLayerRadius        = new double[18]{53.0, 54.6, 56.2, 57.8, 59.4, 61.0, 62.6, 64.2, 65.8, 67.4, 69.0, 70.6, 72.2, 73.8, 75.4, 77.0, 78.6, 80.2};
}


IHoughTransform::~IHoughTransform(){;}

int IHoughTransform::Init()
{
  COMETNamedInfo("IHoughTransform", "//----------------------------------------------------------//");
  COMETNamedInfo("IHoughTransform", "//   Initialize IHoughTransform                               ");
  COMETNamedInfo("IHoughTransform", "//----------------------------------------------------------//");
  return 1;
}

void IHoughTransform::SetHoughTransformVariables(int nIter, int nBins, double nPt, double rhoMax, double rhoMin, int bandWidth, double rad_uncertainty, double refX, double refY){
  fnIter          = nIter;
  fnBins          = nBins;
  fnPt            = nPt;
  fRhoMax         = rhoMax;
  fRhoMin         = rhoMin;
  fBandWidth      = bandWidth;
  fRadUncertainty = rad_uncertainty;
  fRef            = std::make_pair(refX,refY);
}

void IHoughTransform::GetLocalCoordinate(){  
  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
    TVector3 wireend0(fWireEnd0X[i_hit],fWireEnd0Y[i_hit],fWireEnd0Z[i_hit]);	
    TVector3 wireend1(fWireEnd1X[i_hit],fWireEnd1Y[i_hit],fWireEnd1Z[i_hit]);
    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend0, wireend0)){
      continue;
      std::cout << "MisIdentifided wire is detected" << std::endl;}
    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend1, wireend1)){
      continue;
      std::cout << "MisIdentifided wire is detected" << std::endl;}
    fWireEnd0X[i_hit]=wireend1(0);
    fWireEnd0Y[i_hit]=wireend1(1);
    fWireEnd0Z[i_hit]=wireend1(2);
  }  
}

void IHoughTransform::SortHits(){
  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
    // Separate hits in odd and even layer  
    if ((fWireLayerId[i_hit]+1)%2==0){
      fWireEnd0X_even[fnEvenhits]=fWireEnd0X[i_hit];
      fWireEnd0Y_even[fnEvenhits]=fWireEnd0Y[i_hit];
      fnEvenhits++;      
    }
    else if((fWireLayerId[i_hit]+1)%2==1){
      fWireEnd0X_odd[fnOddhits]=fWireEnd0X[i_hit];
      fWireEnd0Y_odd[fnOddhits]=fWireEnd0Y[i_hit];
      fnOddhits++;
    } 

    // Sort layer of hits
    Int_t layerid=fWireLayerId[i_hit];
    if(std::find(fWireIdsPerLayer[layerid].begin(), fWireIdsPerLayer[layerid].end(), fWireId[i_hit]) != fWireIdsPerLayer[layerid].end()){
      continue;
    }
    fWireIdsPerLayer[layerid].push_back(fWireId[i_hit]);
  }
}

void IHoughTransform::ConfigureHitClusters(){
  for (Int_t i_layer=0; i_layer<fNumOfLayers; i_layer++){
    std::sort (fWireIdsPerLayer[i_layer].begin(), fWireIdsPerLayer[i_layer].end(), IntComparison);
    for (Int_t i_ele=0; i_ele<fWireIdsPerLayer[i_layer].size(); i_ele++){
    }
    MakeCluster(fWireIdsPerLayer[i_layer], fTmpClusterSet, fNumOfWiresPerLayer[i_layer]);
  }
  
  // Reduce size 1 Cluster
  for (Int_t i_clu=0; i_clu<fTmpClusterSet.size(); i_clu++){
    if (fTmpClusterSet[i_clu].size()>1){
      fClusterSet.push_back(fTmpClusterSet[i_clu]);
    }
  }  
}

std::vector<std::pair<double,double> > IHoughTransform::MakeOrigins(double rad_uncertainty, std::pair<double,double> ref){
  std::vector<std::pair<double,double> > origins;
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

bool IHoughTransform::PreTrackCut(){
  
  Bool_t ifThereisOddhit=0;
  Bool_t ifThereisEvenhit=0;
  if (fnCALCDCHit<2) return 0; // Require At least more than 2 hits


  //std::cout << fnCALCDCHit << std::endl;
  for (Int_t i_layer=0; i_layer<fnCALCDCHit; i_layer++){      
    if ((fWireLayerId[i_layer]+1)%2==0) {ifThereisEvenhit=1;}
    else if ((fWireLayerId[i_layer]+1)%2==1) {ifThereisOddhit=1;}
  }
  
  if (ifThereisEvenhit==0 || ifThereisOddhit==0) return 0; // Require At least one hit at odd and even layer
  return 1;
}

void IHoughTransform::Process(){

  GetLocalCoordinate();
  SortHits();
  ConfigureHitClusters();

  Int_t deg_index;    
  Int_t rho_index;
  Int_t do_not_use;
  
  Int_t deg_even_index; 
  Int_t rho_even_index;
  Int_t deg_odd_index;
  Int_t rho_odd_index;
  
  Double_t radUncertainty = fRadUncertainty;
  std::pair <Double_t, Double_t> ref = fRef;

  for (Int_t is_even=0; is_even<2; is_even++){    
    for (Int_t it2=0; it2<fnIter; it2++){
      
      std::vector<std::pair<Double_t,Double_t> > origins = MakeOrigins(radUncertainty,ref);     
      std::vector<Int_t> vote_max;
      std::vector<std::pair<Int_t,Int_t> > index_vec; 
      
      for (Int_t it=0; it<origins.size(); it++){
	
	Double_t oriX=origins[it].first; 
	Double_t oriY=origins[it].second;
	
	if  (ifInsideDisk(oriX,oriY)==0){
	  vote_max.push_back(0);
	}	
       
	else if (ifInsideDisk(oriX,oriY)==1){
	  
	  ///////////////////////////////////////
	  ////   Conformal Transformation    ////  
	  ///////////////////////////////////////
	  
	  memset(fConfX,0,sizeof(fConfX));
	  memset(fConfY,0,sizeof(fConfY));	    	  
	  
	  for (Int_t i_hit=0; i_hit<fnCALCDCHit ; i_hit++){   
	    fConfX[i_hit] = ConfTransX(fWireEnd0X[i_hit]-oriX,fWireEnd0Y[i_hit]-oriY);
	    fConfY[i_hit] = ConfTransY(fWireEnd0X[i_hit]-oriX,fWireEnd0Y[i_hit]-oriY);	    
	  }
	  	  
	  ////////////////////////////////////////////
	  ////   Hough Transformation & Voting    ////  
	  ////////////////////////////////////////////
	  
	  Bool_t if_already_vote[fnBins][fnBins][2];
	  std::vector < int > alreadyTransformedId;

	  TH2D *vote = new TH2D("vote_even", "vote_even", fnBins, 0, fnBins, fnBins, 0, fnBins);     
	  
	  //std::cout << "nCALCDCHit " << fnCALCDCHit << std::endl;

	  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
	    memset(if_already_vote,0,sizeof(if_already_vote));

	    if (std::find(alreadyTransformedId.begin(), alreadyTransformedId.end(), fWireId[i_hit]) != alreadyTransformedId.end()) continue;
	    alreadyTransformedId.push_back(fWireId[i_hit]);

	    for (Int_t i_Pt=0 ; i_Pt<fnPt; i_Pt++){
	      
	      Double_t deg = i_Pt/fnPt*180.0;
	      Double_t rho =HoughTrans(fConfX[i_hit],fConfY[i_hit],deg);	  
	      Int_t tmp_deg_index = (deg)*fnBins/180;
	      Int_t tmp_rho_index = (rho+fRhoMax)*fnBins/(fRhoMax-fRhoMin);
	      
	      if (is_even==fWireLayerId[i_hit]%2 && if_already_vote[tmp_deg_index][tmp_rho_index][0]==0){
		if_already_vote[tmp_deg_index][tmp_rho_index][0]=1;
		vote->Fill(tmp_deg_index,tmp_rho_index);

		//std::cout << tmp_deg_index << "   " << tmp_rho_index << std::endl;
		//std::cout << "Filled" << std::endl;
	      }	      
	    }
	  }	  

	  //std::cout << vote->GetMaximumBin() << std::endl;
	  //std::cout << vote->GetBinContent(vote->GetMaximumBin()) << std::endl;
	  //std::cout << vote->GetMaximum() << std::endl;
	  
	  vote_max.push_back(vote->GetMaximum());	      
	  vote->GetMaximumBin(deg_index, rho_index, do_not_use);
	  index_vec.push_back(std::make_pair(deg_index,rho_index));
	  
	  delete vote;	  
	}
      }           
      
   
      Int_t maxIndex=findMaxPoint(vote_max);
      //std::cout << "maxIndex: " <<  maxIndex << std::endl;

      ref.first = origins[maxIndex].first; 
      ref.second = origins[maxIndex].second;
            
      if (is_even==1){	  
	fRef_even.first = ref.first; 
	fRef_even.second = ref.second;
	deg_even_index = index_vec[maxIndex].first; 
	rho_even_index = index_vec[maxIndex].second;
      }
      
      else if (is_even==0){
	fRef_odd.first = ref.first; 
	fRef_odd.second = ref.second;
	deg_odd_index = index_vec[maxIndex].first; 
	rho_odd_index = index_vec[maxIndex].second;
      }      
      radUncertainty = radUncertainty/2.0;      
    }
  }

  Double_t deg_even_eval = 180.0/fnBins*(2*deg_even_index+1)/2;
  Double_t rho_even_eval = fRhoMin+(fRhoMax-fRhoMin)/fnBins*(2*rho_even_index-1)/2;   
  fCx_even = TMath::Cos(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
  fCy_even = TMath::Sin(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
  
  Double_t deg_odd_eval = 180.0/fnBins*(2*deg_odd_index+1)/2;
  Double_t rho_odd_eval = fRhoMin+(fRhoMax-fRhoMin)/fnBins*(2*rho_odd_index-1)/2;
  fCx_odd = TMath::Cos(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);
  fCy_odd = TMath::Sin(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);
  fRad_even = sqrt(pow(fCx_even,2)+pow(fCy_even,2));
  fRad_odd  = sqrt(pow(fCx_odd,2)+pow(fCy_odd,2));   
  fFitpT    = (fRad_even/0.3356+fRad_odd/0.3356)/2;
  fTruthpT  = sqrt(pow(fGenTrPy,2)+pow(fGenTrPz,2));
}

void IHoughTransform::RecognizeHits(){

  /*-------------------------------------------------
    |                                               |
    |    Count Signal Hit Number in Circle Band     |
    |                                               |
    ------------------------------------------------*/

  fOuterR_even = fRad_even +  fBandWidth/2.0;
  fInnerR_even = fRad_even -  fBandWidth/2.0;
  fOuterR_odd  = fRad_odd  +  fBandWidth/2.0;
  fInnerR_odd  = fRad_odd  -  fBandWidth/2.0;
  
  fRecoMaxWireLayerId=0;
  
  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
    if ((fWireLayerId[i_hit]+1)%2 == 0){ // even
      if (ifInsideBand(fWireEnd0X[i_hit],fWireEnd0Y[i_hit],fCx_even+fRef_even.first,fCy_even+fRef_even.second,fOuterR_even,fInnerR_even)==1){
	fRecoWireEnd0X_even[fnRecoHit_even]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y_even[fnRecoHit_even]=fWireEnd0Y[i_hit];
	fnRecoHit_even++;	   
	
	fRecoWireEnd0X[fnRecoHit]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y[fnRecoHit]=fWireEnd0Y[i_hit];
	fRecoWireEnd0Z[fnRecoHit]=fWireEnd0Z[i_hit];
	fRecoWireEnd1X[fnRecoHit]=fWireEnd1X[i_hit];
	fRecoWireEnd1Y[fnRecoHit]=fWireEnd1Y[i_hit];
	fRecoWireEnd1Z[fnRecoHit]=fWireEnd1Z[i_hit];
	fRecoDriftDist[fnRecoHit]=fDriftDist[i_hit];
	fRecoWireLayerId[fnRecoHit]=fWireLayerId[i_hit];
	fRecoWireId[fnRecoHit]=fWireId[i_hit];
	
	if (fRecoWireLayerId[fnRecoHit]>fRecoMaxWireLayerId){
	  fRecoMaxWireLayerId=fRecoWireLayerId[fnRecoHit];
	}
	fnRecoHit++;
      }
    }
    else if ((fWireLayerId[i_hit]+1)%2 == 1){ // odd
      if (ifInsideBand(fWireEnd0X[i_hit],fWireEnd0Y[i_hit],fCx_odd+fRef_odd.first,fCy_odd+fRef_odd.second,fOuterR_odd,fInnerR_odd)==1){
	fRecoWireEnd0X_odd[fnRecoHit_odd]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y_odd[fnRecoHit_odd]=fWireEnd0Y[i_hit];
	fnRecoHit_odd++;
	
	fRecoWireEnd0X[fnRecoHit]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y[fnRecoHit]=fWireEnd0Y[i_hit];
	fRecoWireEnd0Z[fnRecoHit]=fWireEnd0Z[i_hit];
	fRecoWireEnd1X[fnRecoHit]=fWireEnd1X[i_hit];
	fRecoWireEnd1Y[fnRecoHit]=fWireEnd1Y[i_hit];
	fRecoWireEnd1Z[fnRecoHit]=fWireEnd1Z[i_hit];
	fRecoDriftDist[fnRecoHit]=fDriftDist[i_hit];
	fRecoWireLayerId[fnRecoHit]=fWireLayerId[i_hit];
	fRecoWireId[fnRecoHit]=fWireId[i_hit];
	
	if (fRecoWireLayerId[fnRecoHit]>fRecoMaxWireLayerId){
	  fRecoMaxWireLayerId=fRecoWireLayerId[fnRecoHit];
	}
	fnRecoHit++;
      }
    }		
  }    

  /*-----------------------------------------------------------
    |                                                         |
    |           Count Non-Recognized Hits in Cluster          |
    |                                                         |
    ----------------------------------------------------------*/
  
  std::vector <Int_t> WireIdsToBeReco;
  
  for (Int_t i_clu=0; i_clu<fClusterSet.size(); i_clu++){
    Bool_t isItRecoCluster=0;
    
    for(Int_t i_ele=0; i_ele<fClusterSet[i_clu].size(); i_ele++){	  
      if (ifInsideArray(fClusterSet[i_clu][i_ele], fRecoWireId, fnRecoHit)==1){
	isItRecoCluster=1;
	break;
      }
    }
    
    if (isItRecoCluster==1){
      for(Int_t i_ele=0; i_ele<fClusterSet[i_clu].size(); i_ele++){	  
	if (ifInsideArray(fClusterSet[i_clu][i_ele], fRecoWireId, fnRecoHit)==0){
	  WireIdsToBeReco.push_back(fClusterSet[i_clu][i_ele]);
	}
      }
    }
  }
  
  //std::cout << "Num Of Ids to be Recognized: " << WireIdsToBeReco.size() << std::endl;
  
  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
    if (ifInsideVec(fWireId[i_hit],WireIdsToBeReco)==1){
      
      if ((fWireLayerId[i_hit]+1)%2 == 0){ // even
	fRecoWireEnd0X_even[fnRecoHit_even]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y_even[fnRecoHit_even]=fWireEnd0Y[i_hit];
	fnRecoHit_even++;	   
	
	fRecoWireEnd0X[fnRecoHit]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y[fnRecoHit]=fWireEnd0Y[i_hit];
	fRecoWireEnd0Z[fnRecoHit]=fWireEnd0Z[i_hit];
	fRecoWireEnd1X[fnRecoHit]=fWireEnd1X[i_hit];
	fRecoWireEnd1Y[fnRecoHit]=fWireEnd1Y[i_hit];
	fRecoWireEnd1Z[fnRecoHit]=fWireEnd1Z[i_hit];
	
	fRecoDriftDist[fnRecoHit]=fDriftDist[i_hit];
	fRecoWireLayerId[fnRecoHit]=fWireLayerId[i_hit];
	fRecoWireId[fnRecoHit]=fWireId[i_hit];
	
	if (fRecoWireLayerId[fnRecoHit]>fRecoMaxWireLayerId){
	  fRecoMaxWireLayerId=fRecoWireLayerId[fnRecoHit];
	}
	fnRecoHit++;
	
      }
      else if ((fWireLayerId[i_hit]+1)%2 == 1){ // odd
	fRecoWireEnd0X_odd[fnRecoHit_odd]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y_odd[fnRecoHit_odd]=fWireEnd0Y[i_hit];
	fnRecoHit_odd++;
	
	fRecoWireEnd0X[fnRecoHit]=fWireEnd0X[i_hit];
	fRecoWireEnd0Y[fnRecoHit]=fWireEnd0Y[i_hit];
	fRecoWireEnd0Z[fnRecoHit]=fWireEnd0Z[i_hit];
	fRecoWireEnd1X[fnRecoHit]=fWireEnd1X[i_hit];
	fRecoWireEnd1Y[fnRecoHit]=fWireEnd1Y[i_hit];
	fRecoWireEnd1Z[fnRecoHit]=fWireEnd1Z[i_hit];
	fRecoDriftDist[fnRecoHit]=fDriftDist[i_hit];
	fRecoWireLayerId[fnRecoHit]=fWireLayerId[i_hit];
	fRecoWireId[fnRecoHit]=fWireId[i_hit];
	
	if (fRecoWireLayerId[fnRecoHit]>fRecoMaxWireLayerId){
	  fRecoMaxWireLayerId=fRecoWireLayerId[fnRecoHit];
	}
	fnRecoHit++;	    
      }			  
    } 
  }  

  /*--------------------
    |                  |                                      
    |    CL3 Check     |
    |                  |
    -------------------*/
  
  fAbsCx_even=fCx_even+fRef_even.first;
  fAbsCy_even=fCy_even+fRef_even.second;
  fAbsCx_odd =fCx_odd +fRef_odd.first;
  fAbsCy_odd =fCy_odd +fRef_odd.second;      
  std::vector< Int_t > domain1_layer;
  std::vector< TVector3 > domain1_wireend0;
  std::vector< TVector3 > domain1_wireend1;
  std::vector< Int_t > domain2_layer;
  std::vector< TVector3 > domain2_wireend0;
  std::vector< TVector3 > domain2_wireend1;
  
  Double_t slope_even = fCy_even/fCx_even;
  Double_t slope_odd  = fCy_odd/fCx_odd;
  
  for (Int_t i_reco=0; i_reco<fnRecoHit; i_reco++){
    TVector3 wireend0(fRecoWireEnd0X[i_reco],fRecoWireEnd0Y[i_reco],fRecoWireEnd0Z[i_reco]);
    TVector3 wireend1(fRecoWireEnd1X[i_reco],fRecoWireEnd1Y[i_reco],fRecoWireEnd1Z[i_reco]);
    
    if ((fRecoWireLayerId[i_reco]+1)%2 == 0){ //even
      
      Double_t cross=fRecoWireEnd0X[i_reco]*fAbsCy_even-fRecoWireEnd0Y[i_reco]*fAbsCx_even;
      Double_t LHS  =fRecoWireEnd0Y[i_reco];
      Double_t RHS  =slope_even*(fRecoWireEnd0X[i_reco]-fAbsCx_even)+fAbsCy_even;
      if (cross<=0){
	domain1_layer.push_back(fRecoWireLayerId[i_reco]);
	domain1_wireend0.push_back(wireend0);
	domain1_wireend1.push_back(wireend1);
      }
      else if (cross>0){
	domain2_layer.push_back(fRecoWireLayerId[i_reco]);
	domain2_wireend0.push_back(wireend0);
	domain2_wireend1.push_back(wireend1);
      }
    }
    
    else if ((fRecoWireLayerId[i_reco]+1)%2 == 1){ //odd
      Double_t cross=fRecoWireEnd0X[i_reco]*fAbsCy_odd-fRecoWireEnd0Y[i_reco]*fAbsCx_odd;
      Double_t LHS  =fRecoWireEnd0Y[i_reco];
      Double_t RHS  =slope_odd*(fRecoWireEnd0X[i_reco]-fAbsCx_odd)+fAbsCy_odd;
      if (cross<=0){
	domain1_layer.push_back(fRecoWireLayerId[i_reco]);
	domain1_wireend0.push_back(wireend0);
	domain1_wireend1.push_back(wireend1);
      }
      else if (cross>0){
	domain2_layer.push_back(fRecoWireLayerId[i_reco]);
	domain2_wireend0.push_back(wireend0);
	domain2_wireend1.push_back(wireend1);
      }
    }
  }
  
  if (ifInsideVec(0,domain1_layer)==1 && ifInsideVec(1,domain1_layer)==1 && ifInsideVec(2,domain1_layer)==1){
    if (ifInsideVec(0,domain2_layer)==1 && ifInsideVec(1,domain2_layer)==1 && ifInsideVec(2,domain2_layer)==1){
      fRecoCL3=1;
    }
  }

  /*--------------------------------
    |                              |
    |    Charge Identification     |
    |                              |
    -------------------------------*/
       
  std::vector < Int_t > ScintRecoIndex;
  std::vector < Int_t > CherenRecoIndex;  // Indices where Hough Circles are passing through
  std::vector < Bool_t > isIndexClockwise;
  
  for (int i=0; i<fPairCandidates.size(); i++){
    std::vector < int > ScintIndex  = fPairCandidates[i].first;
    std::vector < int > CherenIndex = fPairCandidates[i].second;
    for (int i_sci=0; i_sci<ScintIndex.size(); i_sci++){
      int index=ScintIndex[i_sci];
      if (index<fCTHSegNum){   // Up 
	fScintUpIndex.push_back(index);
      } 
      if (index>=fCTHSegNum){  // Down
	index=(Int_t(1.5*fCTHSegNum)-index%fCTHSegNum)%fCTHSegNum;
	fScintDownIndex.push_back(index);
      } 

      std::pair<Double_t, Double_t> Scintx1y1 = std::make_pair(-fScintWidth/2.,-fScintHeight/2.);
      std::pair<Double_t, Double_t> Scintx1y2 = std::make_pair(-fScintWidth/2.,fScintHeight/2.);
      std::pair<Double_t, Double_t> Scintx2y1 = std::make_pair(fScintWidth/2.,-fScintHeight/2.);
      std::pair<Double_t, Double_t> Scintx2y2 = std::make_pair(fScintWidth/2.,fScintHeight/2.);
      Double_t ScintRotAngle=-(90-(index+1/2.)*360./fCTHSegNum+fScintTiltAngle);
      Double_t xTrans=fScintRad*TMath::Cos((index+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
      Double_t yTrans=fScintRad*TMath::Sin((index+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
      
      Scintx1y1=getRotatedXY(Scintx1y1, ScintRotAngle);
      Scintx1y2=getRotatedXY(Scintx1y2, ScintRotAngle);
      Scintx2y1=getRotatedXY(Scintx2y1, ScintRotAngle);
      Scintx2y2=getRotatedXY(Scintx2y2, ScintRotAngle);
      
      Scintx1y1=std::make_pair(Scintx1y1.first+xTrans, Scintx1y1.second+yTrans);
      Scintx1y2=std::make_pair(Scintx1y2.first+xTrans, Scintx1y2.second+yTrans);
      Scintx2y1=std::make_pair(Scintx2y1.first+xTrans, Scintx2y1.second+yTrans);
      Scintx2y2=std::make_pair(Scintx2y2.first+xTrans, Scintx2y2.second+yTrans);
      
      std::pair <Double_t, Double_t> ref1 = std::make_pair((Scintx1y1.first+Scintx1y2.first)/2,(Scintx1y1.second+Scintx1y2.second)/2);
      std::pair <Double_t, Double_t> ref2 = std::make_pair((Scintx2y1.first+Scintx2y2.first)/2,(Scintx2y1.second+Scintx2y2.second)/2);
      std::pair <Double_t, Double_t> refMid = std::make_pair((ref1.first+ref2.first)/2,(ref1.second+ref2.second)/2);
      Double_t avg_cX=(fAbsCx_even+fAbsCx_odd)/2;
      Double_t avg_cY=(fAbsCy_even+fAbsCy_odd)/2;
      
      if (ifCircleIsPassing(fRad_even, fAbsCx_even, fAbsCy_even, ref1, ref2)==1 || ifCircleIsPassing(fRad_odd, fAbsCx_odd, fAbsCy_odd, ref1, ref2)==1){
	ScintRecoIndex.push_back(index);
	Double_t cross = refMid.first*avg_cY-avg_cX*refMid.second;
	if (cross>=0){
	  isIndexClockwise.push_back(1);
	}
	else if (cross<0){
	  isIndexClockwise.push_back(0);
	}
      }      
    }
    
    for (int i_che=0; i_che<CherenIndex.size(); i_che++){
      int index=CherenIndex[i_che];
      if (index<fCTHSegNum){  //Up 
	fCherenUpIndex.push_back(index);
      } 
      if (index>=fCTHSegNum){  //Down 
	index=(Int_t(1.5*fCTHSegNum)-index%fCTHSegNum)%fCTHSegNum;
	fCherenDownIndex.push_back(index);
      }

      std::pair<Double_t, Double_t> Cherenx1y1 = std::make_pair(-fCherenWidth/2.,-fCherenHeight/2.);
      std::pair<Double_t, Double_t> Cherenx1y2 = std::make_pair(-fCherenWidth/2.,fCherenHeight/2.);
      std::pair<Double_t, Double_t> Cherenx2y1 = std::make_pair(fCherenWidth/2.,-fCherenHeight/2.);
      std::pair<Double_t, Double_t> Cherenx2y2 = std::make_pair(fCherenWidth/2.,fCherenHeight/2.);
      Double_t CherenRotAngle=-(90-(index+1/2.)*360./fCTHSegNum+fCherenTiltAngle);
      Double_t xTrans=fCherenRad*TMath::Cos((index+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
      Double_t yTrans=fCherenRad*TMath::Sin((index+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
      
      Cherenx1y1=getRotatedXY(Cherenx1y1, CherenRotAngle);
      Cherenx1y2=getRotatedXY(Cherenx1y2, CherenRotAngle);
      Cherenx2y1=getRotatedXY(Cherenx2y1, CherenRotAngle);
      Cherenx2y2=getRotatedXY(Cherenx2y2, CherenRotAngle);
      
      Cherenx1y1=std::make_pair(Cherenx1y1.first+xTrans, Cherenx1y1.second+yTrans);
      Cherenx1y2=std::make_pair(Cherenx1y2.first+xTrans, Cherenx1y2.second+yTrans);
      Cherenx2y1=std::make_pair(Cherenx2y1.first+xTrans, Cherenx2y1.second+yTrans);
      Cherenx2y2=std::make_pair(Cherenx2y2.first+xTrans, Cherenx2y2.second+yTrans);
      
      std::pair <Double_t, Double_t> ref1 = std::make_pair((Cherenx1y1.first+Cherenx1y2.first)/2,(Cherenx1y1.second+Cherenx1y2.second)/2);
      std::pair <Double_t, Double_t> ref2 = std::make_pair((Cherenx2y1.first+Cherenx2y2.first)/2,(Cherenx2y1.second+Cherenx2y2.second)/2);
      std::pair <Double_t, Double_t> refMid = std::make_pair((ref1.first+ref2.first)/2,(ref1.second+ref2.second)/2);
      Double_t avg_cX=(fAbsCx_even+fAbsCx_odd)/2;
      Double_t avg_cY=(fAbsCy_even+fAbsCy_odd)/2;
      
      if (ifCircleIsPassing(fRad_even, fAbsCx_even, fAbsCy_even, ref1, ref2)==1 || ifCircleIsPassing(fRad_odd, fAbsCx_odd, fAbsCy_odd, ref1, ref2)==1){
	CherenRecoIndex.push_back(index);
	Double_t cross = refMid.first*avg_cY-avg_cX*refMid.second;
	if (cross>=0){
	  isIndexClockwise.push_back(1);
	}
	else if (cross<0){
	  isIndexClockwise.push_back(0);
	}
      }            
    }
  }
  
  if (isIndexClockwise.empty()==1){
    fReco2DCharge=0; // Non-Classified
  }
  else{
    if (std::find(isIndexClockwise.begin(), isIndexClockwise.end(), 0) == isIndexClockwise.end()){ //vector only contains 1
      fReco2DCharge=1;
    }
    else if (std::find(isIndexClockwise.begin(), isIndexClockwise.end(), 1) == isIndexClockwise.end()){ //vector only contains 0
      fReco2DCharge=-1;
    }
    else {
      fReco2DCharge=0;
    }
  }       
}

//void IHoughTransform::GetMasterCoordinate(){
//  
//}

void IHoughTransform::PrintMCStatus(){
  ITracking::PrintMCStatus();
}

void IHoughTransform::PrintResults(){
  std::cout << "----- HoughTransform Results -----" << std::endl;
  std::cout << "Maximum Layer Id: " << fRecoMaxWireLayerId << std::endl;
  std::cout << "RecoCL3         : " << fRecoCL3 << std::endl;
  std::cout << "Reco2DCharge    : " << fReco2DCharge << std::endl; 
}

void IHoughTransform::DrawEvent(TCanvas* canvas){

  /*-----------------------------------------
    |                                         |
    |   Wire Hit Drawing with Hough Circle    |
    |                                         |
    ------------------------------------------*/

  canvas->cd();
  canvas->DrawFrame(-90,-90,90,90);

  // Hough Circles (Even, Odd)
  
  // TEllipse *Outer_circle_even = new TEllipse(fAbsCx_even,fAbsCy_even,fOuterR_even,fOuterR_even);
  // TEllipse *Inner_circle_even = new TEllipse(fAbsCx_even,fAbsCy_even,fInnerR_even,fInnerR_even);
  // TEllipse *Outer_circle_odd = new TEllipse(fAbsCx_odd,fAbsCy_odd,fOuterR_odd,fOuterR_odd);
  // TEllipse *Inner_circle_odd = new TEllipse(fAbsCx_odd,fAbsCy_odd,fInnerR_odd,fInnerR_odd);

  // Inner_circle_even->SetLineColor(0);  
  // Inner_circle_even->SetFillColor(10);
  // Inner_circle_even->Draw();

  // Inner_circle_odd->SetLineColor(0);  
  // Inner_circle_odd->SetFillColor(10);
  // Inner_circle_odd->Draw();

  // Outer_circle_even->SetLineColor(33);  
  // Outer_circle_even->SetFillColor(33);
  // Outer_circle_even->SetFillStyle(4050);
  // Outer_circle_even->Draw();

  // Outer_circle_odd->SetLineColor(45);  
  // Outer_circle_odd->SetFillColor(45);
  // Outer_circle_odd->SetFillStyle(4050);
  // Outer_circle_odd->Draw();
  
  // Wire Layers  
  
  for (Int_t i=0; i<18; i++){
    TEllipse *WireCircle = new TEllipse(0,0,fLayerRadius[i],fLayerRadius[i]);
    WireCircle->SetFillColor(0);
    WireCircle->SetFillStyle(4000);
    WireCircle->SetLineColor(40);
    WireCircle->Draw();
  }
  
  // Stopping Targets
  
  TEllipse *Disk = new TEllipse(0,0,fDiskRad,fDiskRad);
  Disk->SetFillColor(40);
  Disk->SetLineColor(40);
  Disk->Draw();  
  
  // Trigger Hodoscopes
  
  for (Int_t i_cth=0; i_cth<fCTHSegNum; i_cth++){
    std::pair<Double_t, Double_t> Scintx1y1 = std::make_pair(-fScintWidth/2.,-fScintHeight/2.);
    std::pair<Double_t, Double_t> Scintx1y2 = std::make_pair(-fScintWidth/2.,fScintHeight/2.);
    std::pair<Double_t, Double_t> Scintx2y1 = std::make_pair(fScintWidth/2.,-fScintHeight/2.);
    std::pair<Double_t, Double_t> Scintx2y2 = std::make_pair(fScintWidth/2.,fScintHeight/2.);
    Double_t ScintRotAngle=-(90-(i_cth+1/2.)*360./fCTHSegNum+fScintTiltAngle);
    Double_t xTrans=fScintRad*TMath::Cos((i_cth+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
    Double_t yTrans=fScintRad*TMath::Sin((i_cth+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
    
    Scintx1y1=getRotatedXY(Scintx1y1, ScintRotAngle);
    Scintx1y2=getRotatedXY(Scintx1y2, ScintRotAngle);
    Scintx2y1=getRotatedXY(Scintx2y1, ScintRotAngle);
    Scintx2y2=getRotatedXY(Scintx2y2, ScintRotAngle);
    
    Scintx1y1=std::make_pair(Scintx1y1.first+xTrans, Scintx1y1.second+yTrans);
    Scintx1y2=std::make_pair(Scintx1y2.first+xTrans, Scintx1y2.second+yTrans);
    Scintx2y1=std::make_pair(Scintx2y1.first+xTrans, Scintx2y1.second+yTrans);
    Scintx2y2=std::make_pair(Scintx2y2.first+xTrans, Scintx2y2.second+yTrans);
    
    Double_t x[5]={Scintx1y1.first,  Scintx1y2.first,  Scintx2y2.first,  Scintx2y1.first,  Scintx1y1.first};
    Double_t y[5]={Scintx1y1.second, Scintx1y2.second, Scintx2y2.second, Scintx2y1.second, Scintx1y1.second};
    
    TPolyLine *ScintBox=new TPolyLine(5,x,y);
    ScintBox->SetFillColor(0);
    ScintBox->SetFillStyle(4000);
    
    if (ifInsideVec(i_cth,fScintUpIndex)==1 || ifInsideVec(i_cth,fScintDownIndex)==1){
      ScintBox->SetFillColor(8);
      ScintBox->SetFillStyle(1001);
    }
    
    ScintBox->Draw("f");
    ScintBox->Draw();
  }
  for (Int_t i_cth=0; i_cth<fCTHSegNum; i_cth++){
    std::pair<Double_t, Double_t> Cherenx1y1 = std::make_pair(-fCherenWidth/2.,-fCherenHeight/2.);
    std::pair<Double_t, Double_t> Cherenx1y2 = std::make_pair(-fCherenWidth/2.,fCherenHeight/2.);
    std::pair<Double_t, Double_t> Cherenx2y1 = std::make_pair(fCherenWidth/2.,-fCherenHeight/2.);
    std::pair<Double_t, Double_t> Cherenx2y2 = std::make_pair(fCherenWidth/2.,fCherenHeight/2.);
    Double_t CherenRotAngle=-(90-(i_cth+1/2.)*360./fCTHSegNum+fCherenTiltAngle);
    Double_t xTrans=fCherenRad*TMath::Cos((i_cth+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
    Double_t yTrans=fCherenRad*TMath::Sin((i_cth+1/2.)*360./fCTHSegNum*TMath::Pi()/180.);
    
    Cherenx1y1=getRotatedXY(Cherenx1y1, CherenRotAngle);
    Cherenx1y2=getRotatedXY(Cherenx1y2, CherenRotAngle);
    Cherenx2y1=getRotatedXY(Cherenx2y1, CherenRotAngle);
    Cherenx2y2=getRotatedXY(Cherenx2y2, CherenRotAngle);
    
    Cherenx1y1=std::make_pair(Cherenx1y1.first+xTrans, Cherenx1y1.second+yTrans);
    Cherenx1y2=std::make_pair(Cherenx1y2.first+xTrans, Cherenx1y2.second+yTrans);
    Cherenx2y1=std::make_pair(Cherenx2y1.first+xTrans, Cherenx2y1.second+yTrans);
    Cherenx2y2=std::make_pair(Cherenx2y2.first+xTrans, Cherenx2y2.second+yTrans);
    
    Double_t x[5]={Cherenx1y1.first,  Cherenx1y2.first,  Cherenx2y2.first,  Cherenx2y1.first,  Cherenx1y1.first};
    Double_t y[5]={Cherenx1y1.second, Cherenx1y2.second, Cherenx2y2.second, Cherenx2y1.second, Cherenx1y1.second};
    
    TPolyLine *CherenBox=new TPolyLine(5,x,y);
    
    CherenBox->SetFillColor(0);
    CherenBox->SetFillStyle(4000);
    
    if (ifInsideVec(i_cth,fCherenUpIndex)==1 || ifInsideVec(i_cth,fCherenDownIndex)==1){
      CherenBox->SetFillColor(9);
      CherenBox->SetFillStyle(1001);
    }
    
    CherenBox->Draw("f");
    CherenBox->Draw();
    
  }
  
  
  // Hough Circles (Even, Odd)
  
  TEllipse *circle_even = new TEllipse(fAbsCx_even,fAbsCy_even,fRad_even,fRad_even);
  TEllipse *circle_odd = new TEllipse(fAbsCx_odd,fAbsCy_odd,fRad_odd,fRad_odd);
  TEllipse *center_even = new TEllipse(fAbsCx_even,fAbsCy_even,0.1,0.1);
  TEllipse *center_odd = new TEllipse(fAbsCx_odd,fAbsCy_odd,0.1,0.1);
  
  circle_even->SetFillColor(0);
  circle_even->SetFillStyle(4000);
  circle_even->SetLineColor(4);
  circle_even->SetLineWidth(0);
  circle_even->Draw();
  
  circle_odd->SetFillColor(0);
  circle_odd->SetFillStyle(4000);
  circle_odd->SetLineColor(2);      
  circle_odd->SetLineWidth(0);
  circle_odd->Draw();

  // Hits
  
  TGraph *grEvenhits = new TGraph(fnEvenhits, fWireEnd0X_even, fWireEnd0Y_even);
  grEvenhits->SetTitle("Evenhits");
  grEvenhits->SetMarkerStyle(20);
  grEvenhits->SetMarkerSize(1);
  grEvenhits->SetMarkerColor(4);  //even - blue
  grEvenhits->Draw("P");
  
  TGraph *grOddhits = new TGraph(fnOddhits, fWireEnd0X_odd, fWireEnd0Y_odd);
  grOddhits->SetTitle("Oddhits");
  grOddhits->SetMarkerStyle(20);
  grOddhits->SetMarkerSize(1);
  grOddhits->SetMarkerColor(2);  //odd - red
  grOddhits->Draw("P");
  
  // Recognized Hits
  
  TGraph *grRecoEvenhits = new TGraph(fnRecoHit_even, fRecoWireEnd0X_even, fRecoWireEnd0Y_even);
  grRecoEvenhits->SetTitle("RecoEvenhits");  
  grRecoEvenhits->SetMarkerStyle(20);
  grRecoEvenhits->SetMarkerSize(1);
  grRecoEvenhits->SetMarkerColor(38); // Reco_even - right blue
  grRecoEvenhits->Draw("P");
  
  TGraph *grRecoOddhits = new TGraph(fnRecoHit_odd, fRecoWireEnd0X_odd, fRecoWireEnd0Y_odd);
  grRecoOddhits->SetTitle("RecoOddhits");
  grRecoOddhits->SetMarkerStyle(20);
  grRecoOddhits->SetMarkerSize(1);
  grRecoOddhits->SetMarkerColor(46); // Reco_odd - right red
  grRecoOddhits->Draw("P");      

}

std::vector<int> IHoughTransform::GetRecoWireId(){
  std::vector<int> wireId;
  for (int i=0; i<fnRecoHit; i++){
    wireId.push_back(fRecoWireId[i]);
  }
  return wireId;
}

std::vector<double> IHoughTransform::GetRecoDriftDist(){
  std::vector<double> driftDist;
  for (int i=0; i<fnRecoHit; i++){
    driftDist.push_back(fRecoDriftDist[i]);
  }
  return driftDist;
}

int IHoughTransform::Finish(){
  return 1;
}

/////////////////////////
//                     //
// Auxillary Functions //
//                     //
/////////////////////////

bool IntComparison (int i,int j) { return (i<j); }

void MakeCluster(std::vector<int> &WireIds, std::vector<std::vector< int > > &ClusterSet, int WireNumberInLayer)
{
  if (!WireIds.empty()){
    std::vector<int> tmpCluster;
    std::vector<int> firstCluster;
    tmpCluster.push_back(WireIds[0]);
    firstCluster.push_back(WireIds[0]);

    for(int i_id = 1; i_id < WireIds.size(); ++i_id){
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


bool ifInsideDisk(double x, double y){
  if (pow(x,2)+pow(y,2)<100){ return 1;}
  else {return 0;}
}

bool ifInsideBand(double x, double y, double cX, double cY, double outerR, double interR){
  if ((pow(x-cX,2)+pow(y-cY,2))<pow(outerR,2) && (pow(x-cX,2)+pow(y-cY,2))>pow(interR,2)){
    return 1;
  }
  return 0;
}

bool ifInsideVec(int element, std::vector<int> vec){
  for (int i=0; i<vec.size(); i++){
    if (vec[i]==element){
      return 1;
    }
  }
  return 0;
}

bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size){
  for (Int_t i=0; i<arr_size; i++){
    if (reco_ids[i]==wireid) return 1;
  }
  return 0;
}

Int_t findMaxPoint(std::vector<Int_t> vec){
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
