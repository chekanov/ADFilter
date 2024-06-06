/************************************************************
*	Calculates PCA jet width
*	For details see http://arxiv.org/abs/1002.3982
************************************************************/

#include<iostream>
#include "LParticle.h"
#include "CParticle.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include<vector>
#include "SystemOfUnits.h"


double Width(LParticle jet){
	double jetWidth=0;
        double jetPtSum=0;
        vector<CParticle> constit = jet.GetConstituents();
	TLorentzVector TJet = jet.GetP();
        
         for (unsigned int j=0; j< constit.size(); j++) {
                CParticle cp = constit.at(j);
                double Pt = cp.calcPt()/GeV;
                if (Pt > 0) {
                double dPhi = TJet.Phi() - cp.calcPhi();
                if (dPhi>TMath::Pi()) dPhi=2*TMath::Pi()-dPhi;
                double dEta = TJet.PseudoRapidity() - cp.calcEta();
                double dR=TMath::Sqrt(dPhi*dPhi + dEta *dEta);
                jetWidth += dR * Pt;
                jetPtSum += Pt;
                }
        }



   	jetWidth = jetWidth / jetPtSum;
        return jetWidth;
}
