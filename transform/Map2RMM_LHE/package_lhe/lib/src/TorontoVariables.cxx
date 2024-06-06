#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TDecompSVD.h"
#include <vector>
#include <string>
#include <map>
#include "SystemOfUnits.h"
#include "LParticle.h"
#include "CParticle.h"



// rotation
TLorentzVector* RotateAxes(TLorentzVector p, double M[3][3]){
  double px_rot=M[0][0]*(p.Px())+M[0][1]*(p.Py())+M[0][2]*(p.Pz());
  double py_rot=M[1][0]*(p.Px())+M[1][1]*(p.Py())+M[1][2]*(p.Pz());
  double pz_rot=M[2][0]*(p.Px())+M[2][1]*(p.Py())+M[2][2]*(p.Pz());
  //cout<<"px: "<<p->Px()<<endl; 
  TLorentzVector* prot =new TLorentzVector();
  prot->SetPx(px_rot); 
  prot->SetPy(py_rot);
  prot->SetPz(pz_rot);
  prot->SetE(p.E());
  //=TLorentzVector(px_rot,py_rot,pz_rot,p->E());
  // cout<<"rotated px: "<<prot->Px()<<endl; 
  return prot;
} 


// matrix
int CalcRotationMatrix(double nvec[3],double rot_mat[3][3]){
  //clear momentum tensor
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      //cout<<"\n["<<i<<"]["<<j<<"]: "<<rot_mat[i][j]<<endl;

      rot_mat[i][j]=0.;
    }
  }


  double mag3=sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1]+ nvec[2]*nvec[2]);
  // cout<<"mag3: "<<mag3<<endl; 
  double mag2=sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1]);
  //cout<<"mag2: "<<mag2<<endl; 

  if(mag3<=0){
    cout<<"rotation axis is null"<<endl;
    return 0;
  }

  double ctheta0=nvec[2]/mag3;
  double stheta0=mag2/mag3;
  double cphi0 = (mag2>0.) ? nvec[0]/mag2:0.;
  double sphi0 = (mag2>0.) ? nvec[1]/mag2:0.;

  rot_mat[0][0]=-ctheta0*cphi0;
  rot_mat[0][1]=-ctheta0*sphi0;
  rot_mat[0][2]=stheta0;
  rot_mat[1][0]=sphi0;
  rot_mat[1][1]=-1.*cphi0;
  rot_mat[1][2]=0.;
  rot_mat[2][0]=stheta0*cphi0;
  rot_mat[2][1]=stheta0*sphi0;
  rot_mat[2][2]=ctheta0;

  // cout<<"Done rotation"<<endl;

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      //cout<<"\n["<<i<<"]["<<j<<"]: "<<rot_mat[i][j]<<endl;
    }
  }

  return 1;
}



// actual calculation of variables
std::map<std::string, double> TorontoVariables(LParticle Ljet) {


  TMatrixD MomentumTensor(2,2);

  vector<CParticle> constit = Ljet.GetConstituents();
  TLorentzVector jet = Ljet.GetP();


  double phi0=jet.Phi();
  double eta0=jet.Eta();
  double m0 = jet.M();
  //Jet eta and phi
  double nref[3];
  nref[0]=(cos(phi0)/cosh(eta0));
  nref[1]=(sin(phi0)/cosh(eta0));
  nref[2]=tanh(eta0);
  //this is the rotation matrix
  double RotationMatrix[3][3];
  int p = CalcRotationMatrix(nref,RotationMatrix );


  double Iw00=0.,Iw01=0.,Iw10=0.,Iw11=0.;
  double P2Sum = 0.;
  double PtSum = 0.;
  double PxSum = 0.;
  double PySum = 0.;
  double PzSum = 0.;
  double ESum  = 0.;
 
  //ANGULARITY
  double sum_a=0.;
  const double a=-2.;
 

  for (int j=0; j< constit.size(); j++) {
    CParticle cp = constit.at(j);
    double pt= cp.calcPt()/GeV;
    double eta =cp.calcEta();
    double phi= cp.calcPhi();
    double cx = cp.px()/GeV;
    double cy = cp.py()/GeV;
    double cz = cp.pz()/GeV;
    double cE = cp.e()/GeV;
    TLorentzVector clus = TLorentzVector(cx,cy,cz,cE);


    //CENTRALITY   
    PtSum += pt;
    PxSum += cx;
    PySum += cy;
    PzSum += cz;
    ESum  += cE;


    //ANGULARITY
    double e_i = cE;
    double theta_i = jet.Angle(clus.Vect());
    double e_theta_i=e_i*pow(sin(theta_i),a)*pow(1-cos(theta_i),1-a);
    sum_a += e_theta_i;


    //PLANAR FLOW
    TLorentzVector* rotclus=RotateAxes(clus,RotationMatrix);

    double a=1./(e_i*m0);

    MomentumTensor(0,0) += a*rotclus->Px()*rotclus->Px();
    MomentumTensor(0,1) += a*rotclus->Px()*rotclus->Py();
    MomentumTensor(1,0) += a*rotclus->Py()*rotclus->Px();
    MomentumTensor(1,1) += a*rotclus->Py()*rotclus->Py();

    Iw00+= a*rotclus->Px()*rotclus->Px();
    Iw11+= a*rotclus->Py()*rotclus->Py();
    Iw01+= a*rotclus->Px()*rotclus->Py();
    Iw10+= a*rotclus->Py()*rotclus->Px();

  }//consits
  

    //get the jet mass from consitutents
  double jcmass2 = ESum*ESum-PxSum*PxSum-PySum*PySum-PzSum*PzSum;
  double m=-1.;
  if (jcmass2>0.) m= sqrt(jcmass2);
  double jcmass= m;
  //cout<<"Mass from constituents: "<<jcmass<<endl;


  double Centrality = -1;
  if(ESum > 0) Centrality = PtSum/ESum;


  //get eigenvalues of momentum tensor
  TDecompSVD * aSVD = new TDecompSVD(MomentumTensor);
  TVectorD Lambda = aSVD->GetSig();//diagonal


  //try it maually as a test

  double det=Iw00*Iw11-Iw01*Iw10;
  double trace=Iw00+Iw11;
  double pf=(4*det)/(trace*trace);
  //cout<<"Manual planar flow: "<<pf<<endl;

  double Sphericity = -1;
  Sphericity = 2.*Lambda[1]/(Lambda[0]+Lambda[1]);

  double PlanarFlow = -1;
  PlanarFlow = (4*Lambda[0]*Lambda[1])/((Lambda[0]+Lambda[1])*(Lambda[0]+Lambda[1]));

  //cout<<"Planar flow: "<<PlanarFlow<<endl;

  double Angularity = -1;
  Angularity = sum_a/jcmass;


  delete aSVD;
  
  
  std::map<std::string, double> Variables;
 
  Variables["Sphericity"] = Sphericity;
  Variables["Centrality"] = Centrality;
  Variables["PlanarFlow"] = PlanarFlow;
  Variables["Angularity"] = Angularity;
   
  return Variables;
  
}
