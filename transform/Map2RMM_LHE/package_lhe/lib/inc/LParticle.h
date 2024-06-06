#ifndef LParticle_H 
#define LParticle_H

using namespace std;
#include "TMath.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>
#include <map>
#include <string>
#include "CParticle.h"
 

class LParticle: public TObject {
private:
  Int_t m_type;                 // particle type number
  Int_t m_status;               // status code
  Int_t m_charge;             // Particle charge
  Int_t m_parent;             // mother
  TLorentzVector momentum;    // Initial momentum (x,y,z,tot)
  std::vector<double>  parameters;
  std::vector<CParticle>  constituents;

public:


  LParticle();
  LParticle(Double_t px, Double_t py, Double_t pz, Double_t e, Int_t charge);
  LParticle(LParticle* p);

  LParticle(Int_t charge);
  ~LParticle();

  Int_t GetType() {return m_type;};
  Int_t GetStatus(){return m_status;};
  Int_t GetParent(){return m_parent;};
  Int_t GetCharge(){return m_charge;};
  TLorentzVector GetP(){return momentum;};

  void     SetCharge(Int_t q) {m_charge = q;};
  void     SetParent(Int_t q) {m_parent = q;};
  void     SetType(Int_t q) {m_type = q;};
  void     SetStatus(Int_t q) {m_status = q;};
  void     SetP(const TLorentzVector& mom){ momentum = mom; };
  void     SetParameter(double q) { parameters.push_back( q ); };
  void     SetConstituent(CParticle q) {  constituents.push_back( q ); };

  std::vector<double>   GetParameters() { return parameters; };
  std::vector<CParticle>   GetConstituents() { return constituents; };


ClassDef(LParticle,1)// Monte Carlo particle object
};

#endif

