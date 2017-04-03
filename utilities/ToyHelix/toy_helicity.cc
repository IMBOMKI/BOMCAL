#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitter.h>
#include <DAF.h>

#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include <WirePointMeasurement.h>
#include <PlanarMeasurement.h>
#include <SpacepointMeasurement.h>
#include <ProlateSpacepointMeasurement.h>


#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>


#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>


std::string ZeroPadNumber(int num);

int main(int argc, char *argv[]) {


  //Double_t sigma[3]={0.02,0.02,0.3};
  Double_t sigma[3]={0.02,0.02,0.3};

  Int_t fMinIterations=1;
  Int_t fMaxIterations=10;
  Double_t fMinNDF=2;

  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  // init event display
  //genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  // init fitterc
  //genfit::AbsKalmanFitter* kalman = new genfit::KalmanFitterRefTrack(fMaxIterations);
  genfit::AbsKalmanFitter* kalman = new genfit::KalmanFitter(fMaxIterations);
  //genfit::AbsKalmanFitter* kalman=new genfit::DAF();

  TString path = "/sps/hep/comet/beomki/package/GenFit/test/";
  TString geofile = "genfitVac.root";
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import(path+geofile);
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., 10.)); // 15 kGauss  
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  
  std::string padNum=ZeroPadNumber(atoi(argv[1]));
  std::string inName=padNum+".root";

  std::string inDir;
  std::string outDir;
  if (atoi(argv[2])==1){    
    inDir="/sps/hep/comet/beomki/workspace/project/toy_helix/Al27/";
    outDir="/sps/hep/comet/beomki/workspace/project/toy_helix/Al27/Fit/";
  }

  else if (atoi(argv[2])==0){
    inDir="/sps/hep/comet/beomki/workspace/project/toy_helix/Ti48/";
    outDir="/sps/hep/comet/beomki/workspace/project/toy_helix/Ti48/Fit/";
  }

  std::string outName="fit_"+inName;



  TFile *f = TFile::Open(TString(inDir+inName));
  TTree *t = (TTree*)f->Get("trdata");

  TFile *f_fit = TFile::Open(TString(outDir+outName),"recreate");
  TTree *t_fit = new TTree("trdata","tree data");

  
  Int_t nDet;
  Int_t evtId;
  Double_t x_det[100000];
  Double_t y_det[100000];
  Double_t z_det[100000];
  Double_t ini_x;
  Double_t ini_y;
  Double_t ini_z;
  Double_t ini_px;
  Double_t ini_py;
  Double_t ini_pz;
  Double_t iniP;
  Double_t smearedP;
  
  t->SetBranchAddress("nDet", &nDet);
  t->SetBranchAddress("evtId", &evtId);
  t->SetBranchAddress("x_det", x_det);
  t->SetBranchAddress("y_det", y_det);
  t->SetBranchAddress("z_det", z_det);
  t->SetBranchAddress("ini_x", &ini_x);
  t->SetBranchAddress("ini_y", &ini_y);
  t->SetBranchAddress("ini_z", &ini_z);
  t->SetBranchAddress("ini_px", &ini_px);
  t->SetBranchAddress("ini_py", &ini_py);
  t->SetBranchAddress("ini_pz", &ini_pz);
  t->SetBranchAddress("iniP", &iniP);
  t->SetBranchAddress("smearedP",&smearedP);

  Double_t pVal;
  Double_t p_fit;
  Double_t x_fit[100000];
  Double_t y_fit[100000];
  Double_t z_fit[100000];
  Int_t Reco_charge;

  t_fit->Branch("pVal",&pVal,"pVal/D");
  t_fit->Branch("x_fit",x_fit,"x_fit/D");
  t_fit->Branch("y_fit",y_fit,"y_fit/D");
  t_fit->Branch("z_fit",z_fit,"z_fit/D");
  t_fit->Branch("p_fit",&p_fit, "p_fit/D");
  t_fit->Branch("Reco_charge",&Reco_charge,"Reco_charge/I");
  t_fit->Branch("iniP", &iniP, "iniP/D");
  t_fit->Branch("smearedP",&smearedP,"smearedP/D");

  for (Int_t i_evt=0; i_evt<t->GetEntries(); i_evt++){
    Int_t nVirtualPlanes=0;
    Int_t nTotalHits=0;

    t->GetEntry(i_evt);

    genfit::Track *fitTrack(NULL);
    genfit::AbsTrackRep *rep = new genfit::RKTrackRep(11);
    
    TVector3 posInit=TVector3(ini_x*100,ini_y*100,ini_z*100); // m->cm
    TVector3 momInit=TVector3(ini_px, ini_py, ini_pz); //

    
    genfit::MeasuredStateOnPlane stateInit(rep);
    TMatrixDSym covMInit(6);

    /// Position resolutions
    double resPos[3] = {0,0,0};
    resPos[0] = resPos[1] = resPos[2] = 0.01;  /// 1 mm
    for (int ii=0; ii<3; ii++) { covMInit(ii, ii) = resPos[ii]*resPos[ii]; }
    //bool localCoord=false;

    //set covariant matrix from Momentum and angle resolution
    double sigma_p0     = 0.01; ///   1 MeV
    double sigma_angle0 = 0.1;   /// 100 mrad
    TVector3 dirP    = momInit.Unit();
    TVector3 zaxis   = TVector3(0,0,1);
    TVector3 v_axis  = dirP.Cross(zaxis).Unit();
    TVector3 u_axis  = dirP.Cross(v_axis).Unit();
    TVector3 axis[3] = {dirP, v_axis, u_axis};
    TMatrixD hrot(3, 3);
    for(int ii=0; ii<3; ii++) for(int jj=0; jj<3; jj++) hrot(ii,jj) = axis[jj](ii);
    TMatrixD cov0(3, 3);
    cov0(0,0) = pow(sigma_p0, 2.);
    cov0(1,1) = pow(sigma_angle0*momInit.Mag(), 2.);
    cov0(2,2) = pow(sigma_angle0*momInit.Mag(), 2.);
    TMatrixD covrot(3, 3);
    covrot    = hrot*cov0;
    covrot   *= hrot.T();
    for(int ii=0; ii<3; ii++) for(int jj=0; jj<3; jj++) covMInit[ii+3][jj+3] = covrot(ii, jj);
    rep->setPosMomCov(stateInit, posInit, momInit, covMInit);
    /// Set Track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    rep->get6DStateCov(stateInit, seedState, seedCov);
    fitTrack = new genfit::Track(rep, seedState, seedCov);
    std::vector<genfit::TrackPoint*>     trackPoints; // each track point
    std::vector<genfit::AbsMeasurement*> mesHits;     // each measured hit
    
    for (Int_t i_det=0; i_det<20; i_det++){
      
      trackPoints.push_back(new genfit::TrackPoint());
      TVectorD pos(3);
      TMatrixDSym posMatrix(3);
      pos[0]=x_det[i_det]*100;
      pos[1]=y_det[i_det]*100;
      pos[2]=z_det[i_det]*100;
      

      //smearing
      
      for (Int_t i_xyz=0; i_xyz<3; i_xyz++){	
	pos[i_xyz] += gRandom->Gaus(0, sigma[i_xyz]/TMath::Sqrt(3.));      
      }     
      
      for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
          if (row==col) posMatrix(row,col) = std::pow(sigma[row]/TMath::Sqrt(3.),2);
          else posMatrix(row,col) = 0;
        }
      }
      
      mesHits.push_back( new genfit::SpacepointMeasurement(pos, posMatrix, i_det, i_det, trackPoints.at(i_det)) );

      trackPoints.at(i_det)->addRawMeasurement(mesHits.at(i_det));
      fitTrack->insertPoint(trackPoints.at(i_det), nTotalHits++);
      //nHits = mesHits.size();
    }
    
    /// Do the fitting
    kalman->setMinIterations(fMinIterations);
    
    try{
      kalman->processTrack(fitTrack);
    }    
    catch(genfit::Exception& e){e.what(); continue;}
    
    if (!fitTrack->getFitStatus(rep)->isFitConverged()) {
      kalman->processTrack(fitTrack);
    }
    
    if (!fitTrack->getFitStatus(rep)->isFitted()||!fitTrack->getFitStatus(rep)->isFitConverged()) {
      //COMETNamedInfo("IGenFitter", "Fitting is failed...");
      //delete kalman;
      //delete fitTrack;
      //return NULL;
    }
    if (fitTrack->getFitStatus(rep)->getChi2()<=0 ||
	(fitTrack->getFitStatus(rep)->getNdf()-2*nVirtualPlanes)<fMinNDF) {
      //COMETNamedInfo("IGenFitter", "Fit result might be wrong... (chi2,ndf) = (" << 
      //		     fitTrack->getFitStatus(rep)->getChi2() << "," << 
      //	     (fitTrack->getFitStatus(rep)->getNdf()-2*nVirtualPlanes) << ")");
      //delete kalman;
      //delete fitTrack;
      //return NULL;
    }


    genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, rep);
    genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
    const TVectorD& state = kfsop.getState();
    const TMatrixDSym& cov = kfsop.getCov();

    if (state[0]<0) Reco_charge=-1;
    if (state[0]>0) Reco_charge=1;

    p_fit = Reco_charge/state[0];

    pVal = kalman->getPVal(fitTrack,rep);
    t_fit->Fill();
    //display->addEvent(fitTrack);
  
    std::cout << i_evt << std::endl;

    //delete fitTrack;

}// end loop over events
  f_fit->cd();
  f_fit->Write();
  f_fit->Close();
  f->Close();
  delete kalman;
  //display->open();

}


std::string ZeroPadNumber(int num)
{  
  std::stringstream ss;
  
  // the number is converted to string with the help of stringstream
  ss << num; 
  std::string ret;
  ss >> ret;
   
  // Append zero chars
  int str_length = ret.length();
  for (int i = 0; i < 5 - str_length; i++)
    ret = "0" + ret;
  return ret;
  
}
