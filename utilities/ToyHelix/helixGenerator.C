#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"

Double_t xdot(Double_t C, Double_t yi, Double_t y0, Double_t xdot0);
Double_t ydot(Double_t C, Double_t xi, Double_t x0, Double_t ydot0);

Double_t Getk0(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h);
Double_t Getl0(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h);
Double_t Getk1(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h, Double_t k0);
Double_t Getl1(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h);
Double_t Getk2(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h, Double_t k1);
Double_t Getl2(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h);
Double_t Getk3(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h, Double_t k2);
Double_t Getl3(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h);
Double_t GetX(Double_t xi, Double_t k0, Double_t k1, Double_t k2, Double_t k3);
Double_t GetY(Double_t yi, Double_t l0, Double_t l1, Double_t l2, Double_t l3);
Double_t GetDist(Double_t *DIOenergy, Double_t *par);

void helixGenerator(bool AlorTi, string fileName){

  Int_t eventNum=100;

  Bool_t fUseAl=AlorTi; // 0 is for Ti
  Double_t fMassMu=105.66;
  Double_t fMaxDIOE=105;
  Double_t fMinDIOE;
  Double_t fBindingE;
  Double_t fMassN;
  Double_t a5,a6,a7,a8;
  string outputDir;
  if (fUseAl==1){   // Al27
    fMinDIOE=90.29; // 92.29 MeV - 2 MeV
    a5=8.6434E-17;  
    a6=1.16874E-17;
    a7=-1.87828E-19;
    a8=9.16327E-20;
    fBindingE=0.4754;
    fMassN=23273.122;
    //outputDir=/sps/hep/comet/beomki/workspace/project/toy_helix/Al27/;
    outputDir="/home/bomki/ICEDUST/BOMKI_analysis/v999/utilities/ToyHelix/";
  }
  if (fUseAl==0){   // Ti48
    fMinDIOE=96.88; // 98.88 MeV - 2 MeV
    a5=4.44278E-16;
    a6=9.06648E-17;
    a7=-4.26245E-18;
    a8=8.193E-19;
    fBindingE=1.3616;
    fMassN=44646.861;
    //outputDir=/sps/hep/comet/beomki/workspace/project/toy_helix/;
    outputDir="/home/bomki/ICEDUST/BOMKI_analysis/v999/utilities/ToyHelix/";
  }


  TF1 *DIO_dist = new TF1("delta","GetDist",fMinDIOE,fMaxDIOE,10000);
  DIO_dist->SetParameters(fMassMu,fBindingE,fMassN,a5,a6,a7,a8);
  
  Double_t B=1;  // 1T
  Double_t mass=0.000511; // electron mass 0.511 MeV in GeV
  Int_t charge=-1; // electron
  Double_t coef = charge*B/(mass*3.3356);
  Double_t iniE=0.09229; // signal energy 92.29 MeV in GeV
  //Double_t iniP=TMath::Sqrt(iniE*iniE-mass*mass); // Mag of mom in GeV
  Double_t iniP;
  Double_t step=0.00001; // 0.01mm
  Int_t det_interval = 20;


  ////////////////////////////////////////////////////////////////////////


  Double_t x_pos[100000];
  Double_t y_pos[100000];
  Double_t z_pos[100000];
  Double_t x_det[100000];
  Double_t y_det[100000];
  Double_t z_det[100000];
  Double_t ini_px;
  Double_t ini_py;
  Double_t ini_pz;
  Double_t ini_x=0.;
  Double_t ini_y=0.;
  Double_t ini_z=0.;
  Double_t smearedP;

  Int_t nStep=0;
  Int_t nDet=0;
  Int_t evtId=0;
  Int_t MaxIterations=100000;

  TFile *f = TFile::Open(TString(outputDir+fileName), "recreate");
  TTree *t = new TTree("trdata","tree data");
  t->Branch("step", &step, "step/D");
  t->Branch("evtId", &evtId, "evtId/I");
  t->Branch("nStep", &nStep, "nStep/I");
  t->Branch("nDet", &nDet, "nDet/I");
  t->Branch("x_det", x_det, "x_det[nDet]/D");
  t->Branch("y_det", y_det, "y_det[nDet]/D");
  t->Branch("z_det", z_det, "z_det[nDet]/D");
  t->Branch("ini_x", &ini_x, "ini_x/D");
  t->Branch("ini_y", &ini_y, "ini_y/D");
  t->Branch("ini_z", &ini_z, "ini_z/D");
  t->Branch("ini_px", &ini_px, "ini_px/D");
  t->Branch("ini_py", &ini_py, "ini_py/D");
  t->Branch("ini_pz", &ini_pz, "ini_pz/D");
  t->Branch("iniP", &iniP, "iniP/D");
  t->Branch("smearedP", &smearedP, "smearedP/D");

  for (Int_t i_evt=0; i_evt<eventNum; i_evt++){

    iniP = DIO_dist->GetRandom();
    iniP = 0.001 * iniP; // MeV -> GeV
    TVector3 iniPos=TVector3(ini_x,ini_y,ini_z);

    TRandom *r = new TRandom(0);
    Double_t energyLoss = r->Landau(0.0001,0.0001);
    smearedP=iniP-energyLoss;

    while (smearedP<0.08){
      energyLoss = r->Landau(0.0001,0.0001);
      smearedP=iniP-energyLoss;
    }

    TRandom3 *r3 = new TRandom3(0);
    r3->Sphere(ini_px,ini_py,ini_pz,smearedP);    
    Double_t ini_pT = TMath::Sqrt(ini_px*ini_px+ini_py*ini_py);
    


    while (ini_pT>0.95*smearedP || ini_pz>0.95*smearedP){
      r->Sphere(ini_px,ini_py,ini_pz,smearedP);
      ini_pT = TMath::Sqrt(ini_px*ini_px+ini_py*ini_py);
    }


    Double_t xdot0=ini_px/mass;
    Double_t ydot0=ini_py/mass;  
    Double_t zdot0=ini_pz/mass;
    Double_t x0=iniPos(0);
    Double_t y0=iniPos(1);
    Double_t z0=iniPos(2);
    nStep=0;
    nDet=0;
    
    for (Int_t i=0; i<MaxIterations; i++){
      
      if (i==0){
	//x_vel[0]=xdot0;
	//y_vel[0]=ydot0;
	//z_vel[0]=zdot0;
	x_pos[0]=x0;
	y_pos[0]=y0;
	z_pos[0]=z0;
	x_det[0]=x0;
	y_det[0]=y0;
	z_det[0]=z0;
	continue;
      }
      
      //x_vel[i]=xdot(coef, y_pos[i-1],y0,xdot0);
      //y_vel[i]=ydot(coef, x_pos[i-1],x0,ydot0);
      //z_vel[i]=zdot0;
      
      Double_t k0 = Getk0(coef,y_pos[i-1],y0,xdot0,step);
      Double_t l0 = Getl0(coef,x_pos[i-1],x0,ydot0,step);
      Double_t k1 = Getk1(coef,y_pos[i-1],y0,xdot0,step,k0);
      Double_t l1 = Getl1(coef,x_pos[i-1],x0,ydot0,step);
      Double_t k2 = Getk2(coef,y_pos[i-1],y0,xdot0,step,k1);
      Double_t l2 = Getl2(coef,x_pos[i-1],x0,ydot0,step);
      Double_t k3 = Getk3(coef,y_pos[i-1],y0,xdot0,step,k2);
      Double_t l3 = Getl3(coef,x_pos[i-1],x0,ydot0,step);
      x_pos[i]=GetX(x_pos[i-1],k0,k1,k2,k3);
      y_pos[i]=GetY(y_pos[i-1],l0,l1,l2,l3);
      z_pos[i]=z_pos[i-1]+step*zdot0;
      
      nStep++;

      //std::cout << det_meas << std::endl;

      if (i%det_interval==0){
	x_det[i/det_interval]=x_pos[i];
	y_det[i/det_interval]=y_pos[i];
	z_det[i/det_interval]=z_pos[i];
	nDet++;
      }







      
      if (z_pos[i]>5 || z_pos[i]<-5) break;
    }
      
      
    t->Fill();
    evtId++;
    std::cout << evtId << std::endl;
  }
  
  f->Write();
  f->Close();
  

}





Double_t xdot(Double_t C, Double_t yi, Double_t y0, Double_t xdot0){
  return C*(yi-y0)+xdot0;
}

Double_t ydot(Double_t C, Double_t xi, Double_t x0, Double_t ydot0){
  return -C*(xi-x0)+ydot0;
}

Double_t Getk0(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h){
  return h*xdot(C,yi,y0,xdot0);
}

Double_t Getl0(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h){
  return h*ydot(C,xi,x0,ydot0);
}

Double_t Getk1(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h, Double_t k0){
  return h*xdot(C,yi+1./2*k0,y0,xdot0);
}

Double_t Getl1(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h){
  return h*ydot(C,xi+1./2*h,x0,ydot0);
}

Double_t Getk2(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h, Double_t k1){
  return h*xdot(C,yi+1./2*k1,y0,xdot0);
}

Double_t Getl2(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h){
  return h*ydot(C,xi+1./2*h,x0,ydot0);
}

Double_t Getk3(Double_t C, Double_t yi, Double_t y0, Double_t xdot0, Double_t h, Double_t k2){
  return h*xdot(C,yi+k2,y0,xdot0);
}

Double_t Getl3(Double_t C, Double_t xi, Double_t x0, Double_t ydot0, Double_t h){
  return h*ydot(C,xi+h,x0,ydot0);
}

Double_t GetX(Double_t xi,Double_t k0, Double_t k1, Double_t k2, Double_t k3){
  return xi+1./6*(k0+2*k1+2*k2+k3);
}

Double_t GetY(Double_t yi,Double_t l0, Double_t l1, Double_t l2, Double_t l3){
  return yi+1./6*(l0+2*l1+2*l2+l3);
}

Double_t GetDist(Double_t *DIOenergy, Double_t *par){
  Double_t delta = par[0]-par[1]-DIOenergy[0]-TMath::Power(DIOenergy[0],2)/(2*par[2]);
  Double_t a5=par[3];
  Double_t a6=par[4];
  Double_t a7=par[5];
  Double_t a8=par[6];

  Double_t prob = a5*TMath::Power(delta,5)+a6*TMath::Power(delta,6)+a7*TMath::Power(delta,7)+a8*TMath::Power(delta,8);
  return prob;
}
