/**************************************************
*	Calculates PCA eccentricity
*	For details see http://arxiv.org/abs/1002.3982
**************************************************/

using namespace std; 
#include<iostream>
#include "CParticle.h"
#include "LParticle.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "SystemOfUnits.h"


// see implementation details in
// http://ask.metafilter.com/36213/Best-Fit-Ellipse


double eccentricity(LParticle jet) {

   vector<CParticle> constit = jet.GetConstituents();
   unsigned int num=constit.size();
   TLorentzVector TJet = jet.GetP();
 
   double Dphi[num],Deta[num],E[num];  
   double etaSum = 0.; double phiSum = 0.; double eTot = 0.; 

   for (unsigned int j=0; j< num; j++) {
	CParticle cp = constit.at(j);
        E[j]=cp.e()/GeV;
        Dphi[j] = TJet.Phi() - cp.calcPhi();
        // if (Dphi[j]>TMath::Pi()) Dphi[j]=2*TMath::Pi()-Dphi[j]; 
        if(fabs(Dphi[j]-2.*TMath::Pi())< fabs(Dphi[j])) Dphi[j] -= 2. * TMath::Pi(); 
        if(fabs(Dphi[j]+2.*TMath::Pi())< fabs(Dphi[j])) Dphi[j] += 2. * TMath::Pi(); 
        Deta[j] = TJet.PseudoRapidity() - cp.calcEta();
        etaSum = etaSum + Deta[j] * E[j];
        phiSum = phiSum + Dphi[j] * E[j];
        eTot   = eTot + E[j];
   } 

  etaSum = etaSum/eTot; phiSum = phiSum/eTot;
  for(unsigned int j = 0; j< num; j++) {
    Deta[j] = Deta[j]-etaSum;
    Dphi[j] = Dphi[j]-phiSum;
   }


   double X1=0.; 
   double X2=0; 
   for(unsigned int i = 0; i < num; i++) {
    X1 += 2. * E[i] * Deta[i] * Dphi[i];
    X2 += E[i] * ( Dphi[i] * Dphi[i] - Deta[i] * Deta[i] );
   }


   // variance calculations 
   double Theta = .5 * atan( X1/X2 );
   double sinTheta = TMath::Sin(Theta);
   double cosTheta = TMath::Cos(Theta);
   double Theta2 = Theta + 0.5*TMath::Pi();
   double sinThetaPrime = TMath::Sin(Theta2);
   double cosThetaPrime = TMath::Cos(Theta2);



   double VarX = 0.;
   double VarY = 0.;
   for(unsigned int i = 0; i < num; i++) {
     double X=cosTheta*Deta[i] - sinTheta*Dphi[i];    
     double Y=sinTheta*Deta[i] + cosTheta*Dphi[i]; 
     VarX += E[i]*X*X;
     VarY += E[i]*Y*Y;
   }


// another option (atlas?)
/*
    double VarX = 0.;
    double VarY = 0.;
    for(unsigned int i = 0; i < num; i++) {
    double t1=2.*sinTheta*cosTheta*Deta[i]*Dphi[i];
    double t2=cosTheta*cosTheta*Dphi[i]*Dphi[i];
    double t3=sinTheta*sinTheta*Deta[i]*Deta[i];
    VarX += E[i]*(t1+ t2+ t3);
    t1=2.*sinThetaPrime*cosThetaPrime*Deta[i]*Dphi[i];
    t2=cosThetaPrime*cosThetaPrime*Dphi[i]*Dphi[i];
    t3=sinThetaPrime*sinThetaPrime*Deta[i]*Deta[i];
    VarY += E[i]*(t1+t2+t3);
    }
*/


  double VarianceMax = VarX; 
  double VarianceMin = VarY;
  if(VarianceMax < VarianceMin) {
    VarianceMax = VarY;
    VarianceMin = VarX;
  }

  double ECC=1.0 - (VarianceMin/VarianceMax);
 // cout << "2 nd: VarXmax = " << VarianceMax << "   VarXmin=" << VarianceMin << "  ECC=" << ECC1 << endl;



  return ECC;


} 
