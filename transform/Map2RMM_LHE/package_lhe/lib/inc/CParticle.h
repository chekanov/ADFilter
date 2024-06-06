// class to keep info on constituence
// optimized for small size. Use MeV

#ifndef HEADER_CPARTICLE 
#define HEADER_CPARTICLE


using namespace std;
#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include <string> 
#include <sstream>
#include <vector>


class CParticle  : public TObject {


  protected:

       int  Px,Py,Pz,E,Mass;
       int  Charge;
       int  m_pid;
       int  m_status;
       std::vector<int>  parameters;

  public:

    CParticle();
    CParticle(int pX, int pY, int pZ )
                              {Px=pX; Py=pY; Pz=pZ; Mass=0; Charge=0; 
                               E=0; parameters.clear();  };

    CParticle(int pX, int pY, int pZ, int charge )
                              {Px=pX; Py=pY; Pz=pZ; Charge=charge; E=0; 
                               Mass=0;
                               parameters.clear(); };

    CParticle(int pX, int pY, int pZ, int mass, int charge )
                              {Px=pX; Py=pY; Pz=pZ; Mass=mass; Charge=charge; E=0;
                               parameters.clear(); };

    CParticle(int pX, int pY, int pZ, int mass, int charge, int energy )
                          {Px=pX; Py=pY; Pz=pZ; Mass=mass; Charge=charge; E=energy;
                           parameters.clear();  };


   virtual ~CParticle();
   CParticle(CParticle* p);

   void     setParameter(int q) { parameters.push_back( q ); };
   std::vector<int>   getParameters() { return parameters; };


    // get methods
    int  px(){ return Px;};    
    int  py(){ return Py;};    
    int  pz(){ return Pz;};   
    int  e() { return E;};     
    int    id() { return m_pid;};
    void getPxPyPz(int &px, int &py, int &pz) { px=Px; py=Py; pz=Pz;}; 
    int mass()   { return Mass;};
    int charge() { return Charge;};
    int status() { return m_status;};
 
    // - set values 
    void setID(int c) { m_pid=c; };
    void setCharge(int c) { Charge=c; };
    void setStatus(int c) { m_status=c; };
    void setMass(int  m) { Mass=m; };
    void setE(int e) { E=e; };


    // get invariant mass (assuming that e is defined)
    double invMassE();

     // get invariant mass (assuming that masses are defined)
    double invMassM();


    // calculated values -------------------------------
 
    double  calcPt() { return TMath::Sqrt((double)Px*(double)Px+(double)Py*(double)Py); };
    double  calcPhi() { return TMath::ATan2((double)Py,(double)Px); };
    // from 0 to 2pi
    double  calcPhiWrap();
    // from 0 to 360
    double  calcPhiWrapGrad();

    double  calcP()  { return TMath::Sqrt((double)Px*(double)Px+(double)Py*(double)Py+(double)Pz*(double)Pz); };
    double  calcE()  { return TMath::Sqrt((double)Px*(double)Px+(double)Py*(double)Py+(double)Pz*(double)Pz+(double)Mass*(double)Mass); };
    double  calcEta();

  // copy function
    CParticle copy();

    // additive operation
    CParticle operator+(CParticle p);

   // sort in decreaseing order
    bool operator<( const CParticle& rhs ) const { return (double)Px*(double)Px+(double)Py*(double)Py > rhs.Px*rhs.Px+rhs.Py*rhs.Py; }
   // In your code: vector<myclass> collection;
   // std::sort(collection.begin(), collection.end());


    // print 
    void print();
    string toString();

   ClassDef(CParticle,1);  // Integrate this class into ROOT

};




#endif // HEADER_CPARTICLE

