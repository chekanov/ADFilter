#include "LParticle.h"
#include "TMath.h"
     
ClassImp(LParticle)


LParticle::LParticle() {
  m_type=0; 
  m_charge =0;
  m_status=0;
  m_parent=0;
  momentum  = TLorentzVector(0.,0.,0.,0.);   
  parameters.clear(); 
  constituents.clear();
}


LParticle::~LParticle() {

}

LParticle::LParticle(LParticle* part){
  m_type      = part->GetType();
  m_status    = part->GetStatus();
  m_charge    = part->GetCharge();
  m_parent    = part->GetParent();
  momentum    = part->GetP();
  parameters  = part->GetParameters(); 
  constituents = part->GetConstituents();
}



LParticle::LParticle(Int_t charge) {
  m_charge = charge;
  m_status=0;
  m_type=0;
  momentum  = TLorentzVector(0.,0.,0.,0.);
  parameters.clear(); 
  constituents.clear();
}


LParticle::LParticle(Double_t px, Double_t py, Double_t pz, Double_t e, Int_t charge) {
  m_charge = charge;
  m_status=0;
  m_type=0;
  m_parent=0;
  momentum = TLorentzVector(px,py,pz,e);
  parameters.clear();  
  constituents.clear();
}


