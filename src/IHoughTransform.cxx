#include <IHoughTransform.hxx>

#include "HEPUnits.hxx"
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <assert.h>

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

bool IntComparison (int i,int j);
void MakeCluster(std::vector<int> &WireIds, std::vector<std::vector< int> > &ClusterSet, int WireNumberInLayer);

bool ifInsideDisk(double x, double y);
bool ifInsideBand(double x, double y, double cX, double cY, double outerR, double interR);
bool ifInsideVec(int element, std::vector<int> vec);
bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size);
Int_t findMaxPoint(std::vector<Int_t> vec);
std::pair<Double_t, Double_t> getRotatedXY(std::pair<Double_t, Double_t> xy, Double_t deg);
Bool_t ifCircleIsPassing(Double_t rad, Double_t cX, Double_t cY, std::pair<Double_t, Double_t> ref1, std::pair<Double_t, Double_t> ref2);
double getAngle(double x, double y, double cX, double cY);
void circleResidual(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);
TGraph *POCAgr;
std::pair<TVector3, TVector3> getEnterXYPairs(Double_t x1,Double_t y1,Double_t r1, Double_t r2); // For each of domain1 and domain2

IHoughTransform::IHoughTransform(const char* name = "IHoughTransform", const char* title="houghtransform")
  :ITracking(name,title),
   fnIter(3),
   fnBins(100),
   fnPt(300.0),
   fRhoMax(0.02),
   fRhoMin(-0.02),
   fBandWidth(6),
   fRadUncertainty(5), // <- Disk Radius
   //fRadUncertainty(6),
   fRef(std::make_pair(0,0)),
   fnRecoHit(0),
   fnRecoHit_even(0),
   fnRecoHit_odd(0),
   fRecoCL3(0),
   fRecoCL5(0),
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
  memset(fRecoSideHit,0,sizeof(fRecoSideHit));
  memset(fRecoOuterHit,0,sizeof(fRecoOuterHit));
  memset(fRecoInnerHit,0,sizeof(fRecoInnerHit));

  fHitPairs.clear();
  fPOCAs.clear();
  fnPOCA=0;
  memset(fPOCAx,0,sizeof(fPOCAx));
  memset(fPOCAy,0,sizeof(fPOCAy));
  memset(fPOCAz,0,sizeof(fPOCAz));
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
    TVector3 test;
    TVector3 wireend0(fWireEnd0X[i_hit],fWireEnd0Y[i_hit],fWireEnd0Z[i_hit]);	
    TVector3 wireend1(fWireEnd1X[i_hit],fWireEnd1Y[i_hit],fWireEnd1Z[i_hit]);

    //std::cout << "Master         : " << wireend0(0) << "  " << wireend0(1) << "  " << wireend0(2) << std::endl;

    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend0, wireend0)){
      continue;
      std::cout << "MisIdentifided wire is detected" << std::endl;}
    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend1, wireend1)){
      continue;
      std::cout << "MisIdentifided wire is detected" << std::endl;}
    
    //std::cout << "Master -> Local: " << wireend0(0) << "  " << wireend0(1) << "  " << wireend0(2) << std::endl;
    //if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInGlobalCoordinate(wireend0, test)){
    //  continue;
    //  std::cout << "MisIdentifided wire is detected" << std::endl;}
    //  std::cout << "Local -> Master: " << test(0) << "  " << test(1) << "  " << test(2) << std::endl;

    fWireEnd0X[i_hit]=wireend0(0);
    fWireEnd0Y[i_hit]=wireend0(1);
    fWireEnd0Z[i_hit]=wireend0(2);
    fWireEnd1X[i_hit]=wireend1(0);
    fWireEnd1Y[i_hit]=wireend1(1);
    fWireEnd1Z[i_hit]=wireend1(2);
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
	
	     
	if  (ifInsideDisk(oriX,oriY)==0){  // With Inside Disk option
	  vote_max.push_back(0);
	}	
       
	else if (ifInsideDisk(oriX,oriY)==1){
	
	//if (ifInsideDisk(oriX,oriY)>-1){  // Without Inside Disk Option
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
  fTruthpT  = sqrt(pow(fCDCEnterPy,2)+pow(fCDCEnterPz,2));
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
	fRecoDomain[i_reco]=1;
      }
      else if (cross>0){
	domain2_layer.push_back(fRecoWireLayerId[i_reco]);
	domain2_wireend0.push_back(wireend0);
	domain2_wireend1.push_back(wireend1);
	fRecoDomain[i_reco]=2;
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
	fRecoDomain[i_reco]=1;
      }
      else if (cross>0){
	domain2_layer.push_back(fRecoWireLayerId[i_reco]);
	domain2_wireend0.push_back(wireend0);
	domain2_wireend1.push_back(wireend1);
	fRecoDomain[i_reco]=2;
      }
    }
  }
  
  fRecoCL3=0;
  if (ifInsideVec(0,domain1_layer)==1 && ifInsideVec(1,domain1_layer)==1 && ifInsideVec(2,domain1_layer)==1){
    if (ifInsideVec(0,domain2_layer)==1 && ifInsideVec(1,domain2_layer)==1 && ifInsideVec(2,domain2_layer)==1){
      fRecoCL3=1;
    }
  }

  fRecoCL5=0;
  if (ifInsideVec(0,domain1_layer)==1 && ifInsideVec(1,domain1_layer)==1 && ifInsideVec(2,domain1_layer)==1 && ifInsideVec(3,domain1_layer)==1 && ifInsideVec(4,domain1_layer)==1){
    if (ifInsideVec(0,domain2_layer)==1 && ifInsideVec(1,domain2_layer)==1 && ifInsideVec(2,domain2_layer)==1 && ifInsideVec(3,domain2_layer)==1 && ifInsideVec(4,domain2_layer)==1){
      fRecoCL5=1;
    }
  }

  /*------------------------------
    |                            |
    |    Configure Side Hits     |
    |                            |
    -----------------------------*/

  std::map<int,int> Layer;
  std::map<int,int> WireId;
  std::map<int,int> Domain;
  std::map<int,double> Angle;

  for (Int_t i_reco=0; i_reco<fnRecoHit; i_reco++){
    Layer.insert(std::make_pair(i_reco,fRecoWireLayerId[i_reco]));
    Domain.insert(std::make_pair(i_reco,fRecoDomain[i_reco]));
    WireId.insert(std::make_pair(i_reco,fRecoWireId[i_reco]));
    if ((fRecoWireLayerId[i_reco]+1)%2 == 0){      //even	
      Angle.insert(std::make_pair(i_reco,getAngle(fRecoWireEnd0X[i_reco],fRecoWireEnd0Y[i_reco],fAbsCx_even, fAbsCy_even)) );
    }
    else if ((fRecoWireLayerId[i_reco]+1)%2 == 1){ //odd	
      Angle.insert(std::make_pair(i_reco,getAngle(fRecoWireEnd0X[i_reco],fRecoWireEnd0Y[i_reco],fAbsCx_odd, fAbsCy_odd)) );
    }
  }

  std::vector<int> SideWireId;
  std::vector<int> OuterWireId;
  std::vector<int> InnerWireId;

  // Domain 1
  for (Int_t i_dom=1  ; i_dom<=2; i_dom++){  
    for (Int_t i_layer=0; i_layer<=fRecoMaxWireLayerId; i_layer++){
      Double_t maxAngle=0;
      Double_t minAngle=10000;
      Int_t    outerWireId=-1;
      Int_t    innerWireId=-1;
      for (Int_t i=0; i<Layer.size(); i++){
	if (Layer.at(i)==i_layer && Domain.at(i)==i_dom){
	  if (Angle.at(i)>maxAngle){
	    maxAngle  = Angle.at(i);	  
	    outerWireId = WireId.at(i);
	  }	  
	  if (Angle.at(i)<minAngle){
	    minAngle  = Angle.at(i);	  
	    innerWireId = WireId.at(i);
	  }	  
	}
      }    
      if (outerWireId>0) SideWireId.push_back(outerWireId);
      if (outerWireId>0) OuterWireId.push_back(outerWireId);
      if (innerWireId>0) InnerWireId.push_back(innerWireId);
    }
  }

  for (int i_reco=0; i_reco<fnRecoHit; i_reco++){
    int wireId = fRecoWireId[i_reco];
    if (std::find(SideWireId.begin(), SideWireId.end(), wireId) != SideWireId.end()){
      fRecoSideHit[i_reco]=1;
    }
    if (std::find(OuterWireId.begin(), OuterWireId.end(), wireId) != OuterWireId.end()){
      fRecoOuterHit[i_reco]=1;
    }
    if (std::find(InnerWireId.begin(), InnerWireId.end(), wireId) != InnerWireId.end()){
      fRecoInnerHit[i_reco]=1;
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

void IHoughTransform::AddSideHits(){
  //assert(n<=fRecoMaxWireLayerId);
  for (int domain=1; domain<=2; domain++){
    for (int i=0; i<fRecoMaxWireLayerId; i++){ 
   
      struct SingleHit hit_lo, hit_up;
      struct HitPair hitPair;
      
      int i_lo=-1;
      int i_up=-1;
      
      for (int i_hit=0; i_hit<fnRecoHit; i_hit++){
	if      (fRecoWireLayerId[i_hit]==i   && fRecoDomain[i_hit]==domain && fRecoInnerHit[i_hit]==1) i_lo=i_hit;
	else if (fRecoWireLayerId[i_hit]==i+1 && fRecoDomain[i_hit]==domain && fRecoOuterHit[i_hit]==1) i_up=i_hit;
      }
      
      if (i_lo==-1 || i_up==-1) continue;
      
      TVector3 wireEnd0_lo = TVector3(fRecoWireEnd0X[i_lo], fRecoWireEnd0Y[i_lo], fRecoWireEnd0Z[i_lo]);
      TVector3 wireEnd1_lo = TVector3(fRecoWireEnd1X[i_lo], fRecoWireEnd1Y[i_lo], fRecoWireEnd1Z[i_lo]);
      
      TVector3 wireEnd0_up = TVector3(fRecoWireEnd0X[i_up], fRecoWireEnd0Y[i_up], fRecoWireEnd0Z[i_up]);
      TVector3 wireEnd1_up = TVector3(fRecoWireEnd1X[i_up], fRecoWireEnd1Y[i_up], fRecoWireEnd1Z[i_up]);
      
      hit_lo.wireEnd0  = wireEnd0_lo;
      hit_lo.wireEnd1  = wireEnd1_lo;
      hit_lo.driftDist = fRecoDriftDist[i_lo];
      hit_lo.wireId    = fRecoWireId[i_lo];
      hit_lo.layerId   = fRecoWireLayerId[i_lo];
      hit_lo.domain    = fRecoDomain[i_lo];
      
      hit_up.wireEnd0  = wireEnd0_up;
      hit_up.wireEnd1  = wireEnd1_up;
      hit_up.driftDist = fRecoDriftDist[i_up];
      hit_up.wireId    = fRecoWireId[i_up];
      hit_up.layerId   = fRecoWireLayerId[i_up];
      hit_up.domain    = fRecoDomain[i_up];
      
      //std::cout << wireEnd0_lo(0) << "  " << wireEnd0_lo(1) << "  " << wireEnd0_lo(2) << std::endl;
      //std::cout << wireEnd1_lo(0) << "  " << wireEnd1_lo(1) << "  " << wireEnd1_lo(2) << std::endl;
      //std::cout << wireEnd0_up(0) << "  " << wireEnd0_up(1) << "  " << wireEnd0_up(2) << std::endl;
      //std::cout << wireEnd1_up(0) << "  " << wireEnd1_up(1) << "  " << wireEnd1_up(2) << std::endl;

      TVector3 POCA = GetPOCAofTwoWires(hit_lo.wireEnd0, hit_lo.wireEnd1, hit_up.wireEnd0, hit_up.wireEnd1);
      TVector3 cVec = GetVectorCrossingCenter(hit_lo.wireEnd0, hit_lo.wireEnd1, hit_up.wireEnd0, hit_up.wireEnd1,POCA);
      
      hitPair.h1 = hit_lo;
      hitPair.h2 = hit_up;
      hitPair.cV = cVec;
      
      fHitPairs.push_back(hitPair);
      fPOCAs.push_back(POCA);
      
      //std::cout << "Layer: <" << i << " " << i+1 << ">   Domain: " << domain << std::endl;
      //std::cout << "POCA:  " << POCA(0) << "  " << POCA(1) << "  " << POCA(2) << std::endl;
    }
  }
  
  for (Int_t i=0; i<fPOCAs.size(); i++){    
    fPOCAx[fnPOCA]=fPOCAs.at(i)(0);
    fPOCAy[fnPOCA]=fPOCAs.at(i)(1);
    fPOCAz[fnPOCA]=fPOCAs.at(i)(2);

    //std::cout << fPOCAx[fnPOCA] << "  " << fPOCAy[fnPOCA] << "  " << fPOCAz[fnPOCA] << std::endl;
    fnPOCA++;
  }
}

void IHoughTransform::TuneRadiusWithPOCAs(){
  AddSideHits();
  POCAgr = new TGraph(fnPOCA,fPOCAx,fPOCAy);
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetFCN(circleResidual);
  
  Double_t seedCx   = (fAbsCx_even+fAbsCx_odd)/2;
  Double_t seedCy   = (fAbsCy_odd+fAbsCy_odd)/2;
  Double_t seedRad  = (fRad_even+fRad_odd)/2;

  fitter->SetParameter(0, "seedCx", seedCx,  0.1, seedCx-10,seedCx+10);
  fitter->SetParameter(1, "seedCy", seedCy,  0.1, seedCy-10,seedCy+10);
  fitter->SetParameter(2, "seedRad",seedRad, 0.1, seedRad-6,seedRad+6);

  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);

  fAbsCx_Reseeded = fitter->GetParameter(0);
  fAbsCy_Reseeded = fitter->GetParameter(1);
  fFitR_Reseeded  = fitter->GetParameter(2);
  fFitpT_Reseeded = (fitter->GetParameter(2))/0.3356;
  std::cout << "Seed Value: " << seedCx << "  " << seedCy << "  " << seedRad/0.3356 << std::endl;
  std::cout << "Result    : " << fAbsCx_Reseeded << "  " << fAbsCy_Reseeded << "  " << fFitpT_Reseeded << std::endl;

  fXY_domain1 = getEnterXYPairs(fAbsCx_Reseeded,fAbsCy_Reseeded,fFitR_Reseeded,53).first;
  fXY_domain2 = getEnterXYPairs(fAbsCx_Reseeded,fAbsCy_Reseeded,fFitR_Reseeded,53).second;
  fEnterX_domain1 = fXY_domain1(0);
  fEnterY_domain1 = fXY_domain1(1);
  fEnterX_domain2 = fXY_domain2(0);
  fEnterY_domain2 = fXY_domain2(1);

  TVector3 radVec_domain1 = fXY_domain1-TVector3(fAbsCx_Reseeded,fAbsCy_Reseeded,0);
  TVector3 radVec_domain2 = fXY_domain2-TVector3(fAbsCx_Reseeded,fAbsCy_Reseeded,0);
  radVec_domain1.SetMag(fFitpT_Reseeded);
  radVec_domain2.SetMag(fFitpT_Reseeded);
  fPxPy_domain1 = TVector3(0,0,-1).Cross(radVec_domain1);
  fPxPy_domain2 = TVector3(0,0, 1).Cross(radVec_domain2);

  fEnterPx_domain1 = fPxPy_domain1(0);
  fEnterPy_domain1 = fPxPy_domain1(1);
  fEnterPx_domain2 = fPxPy_domain2(0);
  fEnterPy_domain2 = fPxPy_domain2(1);

  std::cout << "Domain 1 Enter: " << fEnterX_domain1 << "  " << fEnterY_domain1 << std::endl;
  std::cout << "Domain 2 Enter: " << fEnterX_domain2 << "  " << fEnterY_domain2 << std::endl;
}

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
  /*
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
  */


  // Side Hits
  /*
  double fOuterWireEnd0X[1000];
  double fOuterWireEnd0Y[1000];  
  int    fnOuterHit=0;
  for (int i_reco=0; i_reco<fnRecoHit; i_reco++){
    if (fRecoOuterHit[i_reco]==1){
      fOuterWireEnd0X[fnOuterHit]=fRecoWireEnd0X[i_reco];
      fOuterWireEnd0Y[fnOuterHit]=fRecoWireEnd0Y[i_reco];
      fnOuterHit++;
    }
  }

  TGraph *grOuterHits = new TGraph(fnOuterHit, fOuterWireEnd0X, fOuterWireEnd0Y);
  grOuterHits->SetTitle("OuterHits");
  grOuterHits->SetMarkerStyle(20);
  grOuterHits->SetMarkerSize(1);
  grOuterHits->SetMarkerColor(5); // Outer Hit - Yellow 
  grOuterHits->Draw("P");        

  double fInnerWireEnd0X[1000];
  double fInnerWireEnd0Y[1000];  
  int    fnInnerHit=0;
  for (int i_reco=0; i_reco<fnRecoHit; i_reco++){
    if (fRecoInnerHit[i_reco]==1){
      fInnerWireEnd0X[fnInnerHit]=fRecoWireEnd0X[i_reco];
      fInnerWireEnd0Y[fnInnerHit]=fRecoWireEnd0Y[i_reco];
      fnInnerHit++;
    }
  }

  TGraph *grInnerHits = new TGraph(fnInnerHit, fInnerWireEnd0X, fInnerWireEnd0Y);
  grInnerHits->SetTitle("InnerHits");
  grInnerHits->SetMarkerStyle(20);
  grInnerHits->SetMarkerSize(1);
  grInnerHits->SetMarkerColor(8); // Inner Hit - Green
  grInnerHits->Draw("P");        
  */


  // POCAs
  POCAgr->SetTitle("POCAs");
  POCAgr->SetMarkerStyle(20);
  POCAgr->SetMarkerSize(1);
  POCAgr->SetMarkerColor(37); 
  POCAgr->Draw("P");  

  /*
  TGraph *ptEnterXY_domain1 = new TGraph(1, &fEnterX_domain1, &fEnterY_domain1);
  ptEnterXY_domain1->SetTitle("EnterXY_domain1");
  ptEnterXY_domain1->SetMarkerStyle(20);
  ptEnterXY_domain1->SetMarkerSize(1);
  ptEnterXY_domain1->SetMarkerColor(30);
  ptEnterXY_domain1->Draw("P");        

  TGraph *ptEnterXY_domain2 = new TGraph(1, &fEnterX_domain2, &fEnterY_domain2);
  ptEnterXY_domain2->SetTitle("EnterXY_domain2");
  ptEnterXY_domain2->SetMarkerStyle(20);
  ptEnterXY_domain2->SetMarkerSize(1);
  ptEnterXY_domain2->SetMarkerColor(20);
  ptEnterXY_domain2->Draw("P");       

  TGraph *ptCenterXY = new TGraph(1, &fAbsCx_Reseeded, &fAbsCy_Reseeded);
  ptCenterXY->SetTitle("CenterXY");
  ptCenterXY->SetMarkerStyle(20);
  ptCenterXY->SetMarkerSize(1);
  ptCenterXY->SetMarkerColor(1);
  ptCenterXY->Draw("P");       
  */

  TEllipse *Circle_Reseeded = new TEllipse(fAbsCx_Reseeded,fAbsCy_Reseeded,fFitR_Reseeded,fFitR_Reseeded);
  Circle_Reseeded->SetFillColor(0);
  Circle_Reseeded->SetFillStyle(4000);
  Circle_Reseeded->SetLineColor(12);
  Circle_Reseeded->Draw();
  
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

std::vector<int> IHoughTransform::GetRecoDomain(){
  std::vector<int> domain;
  for (int i=0; i<fnRecoHit; i++){
    domain.push_back(fRecoDomain[i]);
  }
  return domain;
}

std::vector<bool> IHoughTransform::GetRecoSideHit(){
  std::vector<bool> side;
  for (int i=0; i<fnRecoHit; i++){
    side.push_back(fRecoSideHit[i]);
  }
  return side;
}

std::vector<bool> IHoughTransform::GetRecoOuterHit(){
  std::vector<bool> side;
  for (int i=0; i<fnRecoHit; i++){
    side.push_back(fRecoOuterHit[i]);
  }
  return side;
}

std::vector<bool> IHoughTransform::GetRecoInnerHit(){
  std::vector<bool> side;
  for (int i=0; i<fnRecoHit; i++){
    side.push_back(fRecoInnerHit[i]);
  }
  return side;
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

double getAngle(double x, double y, double cX, double cY){
  double cross = TMath::Abs(x*cY - y*cX);
  double v1Mag = TMath::Sqrt(x*x+y*y);
  double v2Mag = TMath::Sqrt(cX*cX+cY*cY);
  double angle = TMath::ASin(cross/(v1Mag*v2Mag));
  return angle;
}

void circleResidual(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  //minimisation function computing the sum of squares of residuals
  Int_t np = POCAgr->GetN();
  f = 0;
  Double_t *x = POCAgr->GetX();
  Double_t *y = POCAgr->GetY();
  for (Int_t i=0;i<np;i++) {
    Double_t u = x[i] - par[0];
    Double_t v = y[i] - par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    f += dr*dr;
  }
}

std::pair<TVector3, TVector3> getEnterXYPairs(Double_t x1,Double_t y1,Double_t r1, Double_t r2){
  TVector3 XY_domain1;
  TVector3 XY_domain2;
  Double_t x,y,R,SqTerm,cross;
  R=sqrt(pow(x1,2)+pow(y1,2));
  SqTerm=2*(pow(r1,2)+pow(r2,2))/pow(R,2)-pow((pow(r1,2)-pow(r2,2)),2)/pow(R,4)-1;
  if (SqTerm>0)      SqTerm = sqrt(SqTerm);
  else if (SqTerm<0) return std::make_pair(XY_domain1,XY_domain2); // return Nothing

  x=0.5*x1+0.5*(pow(r1,2)-pow(r2,2))/pow(R,2)*(-x1)+0.5*(-y1)*SqTerm;
  y=0.5*y1+0.5*(pow(r1,2)-pow(r2,2))/pow(R,2)*(-y1)+0.5* (x1)*SqTerm;

  cross = x*y1-y*x1;
  if (cross<=0)     XY_domain1=TVector3(x,y,0.);
  else if (cross>0) XY_domain2=TVector3(x,y,0.);

  x=0.5*x1+0.5*(pow(r1,2)-pow(r2,2))/pow(R,2)*(-x1)-0.5*(-y1)*SqTerm;
  y=0.5*y1+0.5*(pow(r1,2)-pow(r2,2))/pow(R,2)*(-y1)-0.5* (x1)*SqTerm;
  cross = x*y1-y*x1;

  if (cross<=0)     XY_domain1=TVector3(x,y,0.);
  else if (cross>0) XY_domain2=TVector3(x,y,0.);

  return std::make_pair(XY_domain1,XY_domain2);
}
