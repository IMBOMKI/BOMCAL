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


bool IntComparison (int i,int j);

void MakeCluster(std::vector<int> &WireIds, std::vector<std::vector< int> > &ClusterSet, int WireNumberInLayer);

IHoughTransform::IHoughTransform(const char* name = "IHoughTransform", const char* title="houghtransform"):ITracking(name,title)
{;}
IHoughTransform::~IHoughTransform(){;}

int IHoughTransform::Init()
{
  COMETNamedInfo("IHoughTransform", "//----------------------------------------------------------//");
  COMETNamedInfo("IHoughTransform", "//   Initialize IHoughTransform                               ");
  COMETNamedInfo("IHoughTransform", "//----------------------------------------------------------//");
  return 1;
}

void IHoughTransform::GetDSCoordinate(){
  /*
  for (Int_t i_hit=0; i_hit<fnCALCDCHit; i_hit++){
    TVector3 wireend0(fWireEnd0X[i_hit],fWireEnd0Y[i_hit],fWireEnd0Z[i_hit]);	
    TVector3 wireend1(fWireEnd1X[i_hit],fWireEnd1Y[i_hit],fWireEnd1Z[i_hit]);
    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend0, wireend0)){
      continue;
      std::cout << "MisIdentifided wire is detected" << std::endl;}
    if(!COMET::IGeomInfo::DetectorSolenoid().GetDetPositionInDSCoordinate(wireend1, wireend1)){
      continue;
      std::cout << "MisIdentifided wire is detected" << std::endl;}
  }
  */
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
