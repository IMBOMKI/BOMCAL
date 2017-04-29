#include <IHoughTransform.hxx>

#include "HEPUnits.hxx"
#include <vector>
#include <iostream>

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

bool IntComparison (int i,int j);

void MakeCluster(std::vector<int> &WireIds, std::vector<std::vector< int> > &ClusterSet, int WireNumberInLayer);

bool ifInsideDisk(double x, double y);
bool ifInsideBand(double x, double y, double cX, double cY, double outerR, double interR);
bool ifInsideVec(int element, std::vector<int> vec);
Int_t findMaxPoint(std::vector<Int_t> vec);

IHoughTransform::IHoughTransform(const char* name = "IHoughTransform", const char* title="houghtransform")
  :ITracking(name,title),
   fnIter(3),
   fnBins(100),
   fnPt(300.0),
   fRhoMax(0.02),
   fRhoMin(-0.02),
   fBandWidth(8),
   fRadUncertainty(5),
   fRef(std::make_pair(0,0))
{;}
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
  std::pair <Double_t,Double_t> ref_even;
  Int_t deg_odd_index;
  Int_t rho_odd_index;
  std::pair <Double_t,Double_t> ref_odd;
  TH2F *ref_even_dist = new TH2F("ref_even_dist", "ref_even_dist", 10, -10, 10, 10, -10, 10);
  TH2F *ref_odd_dist = new TH2F("ref_odd_dist", "ref_odd_dist", 10, -10, 10, 10, -10, 10);
  
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

	    //std::cout << fWireEnd0X[i_hit] << "   " << fWireEnd0Y[i_hit] << "   " << fWireEnd0Z[i_hit] << std::endl;
	    //std::cout << fConfX[i_hit] << "   " << fConfY[i_hit] << "   "  << std::endl;	   

	  }
	  	  
	  ////////////////////////////////////////////
	  ////   Hough Transformation & Voting    ////  
	  ////////////////////////////////////////////
	  
	  Bool_t if_already_vote[fnBins][fnBins][2];
	  TH2D *vote = new TH2D("vote_even", "vote_even", fnBins, 0, fnBins, fnBins, 0, fnBins);     
	  
	  //std::cout << "nCALCDCHit " << fnCALCDCHit << std::endl;

	  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
	    memset(if_already_vote,0,sizeof(if_already_vote));
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
	ref_even.first = ref.first; 
	ref_even.second = ref.second;
	deg_even_index = index_vec[maxIndex].first; 
	rho_even_index = index_vec[maxIndex].second;
	
	if (it2==fnIter-1){
	  ref_even_dist->Fill(ref_even.first,ref_even.second);}	
      }
      
      else if (is_even==0){
	ref_odd.first = ref.first; 
	ref_odd.second = ref.second;
	deg_odd_index = index_vec[maxIndex].first; 
	rho_odd_index = index_vec[maxIndex].second;
	
	if (it2==fnIter-1){
	  ref_odd_dist->Fill(ref_odd.first,ref_odd.second);}	
      }      
      radUncertainty = radUncertainty/2.0;
      
    }
  }

  std::cout << ref_even.first << "  " << ref_even.second << std::endl;
  std::cout << ref_odd.first  << "  " << ref_odd.second  << std::endl;

  delete ref_even_dist;
  delete ref_odd_dist;
}

void IHoughTransform::GetMasterCoordinate(){
  
}

void IHoughTransform::PrintMCStatus(){
  ITracking::PrintMCStatus();
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
