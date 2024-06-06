// a class to keep a constinuent inside a jet 

#include"CParticle.h"

ClassImp(CParticle) // Integrate this class into ROOT

CParticle::CParticle() {
        Px=0; Py=0; Pz=0; Mass=0; Charge=0;
        parameters.clear();  
        m_pid=0;
        m_status=0;
      } 


CParticle::~CParticle() {
}


// make a copy 
CParticle::CParticle(CParticle* part){
  m_pid    = part->id();
  m_status = part->status();
  Charge  = part->charge();
  Mass    = part->mass();
  Px  = part->px();
  Py  = part->py();
  Pz  = part->pz();
  E  = part->e();
  parameters = part->getParameters();
}



// overloading + for particlses
CParticle CParticle::operator+(CParticle p) {

   CParticle temp;
   temp.Px = Px+p.Px;
   temp.Py = Py+p.Py;
   temp.Pz = Pz+p.Pz;
   temp.E =  E+p.E;
   temp.Charge = Charge+p.Charge;
   temp.Mass   = Mass+p.Mass;
   temp.m_status=m_status;
   temp.m_pid=m_pid;

   return temp;


 }



// copy
CParticle CParticle::copy() {
   CParticle temp;
   temp.Px = Px;
   temp.Py = Py;
   temp.Pz = Pz;
   temp.E =  E;
   temp.Charge = Charge;
   temp.Mass   = Mass;
   temp.m_pid  = m_pid;
   temp.parameters  = parameters;
   temp.m_status  = m_status;
 
   return temp;
  }


// getEta 
double CParticle::calcEta() {
    double  Pt_trk=TMath::Sqrt((double)Px*(double)Px+(double)Py*(double)Py);
    return -log(tan(atan2(Pt_trk,(double)Pz)/2));
}


// getEta
double CParticle::calcPhiWrap() {
   
    double  phi=(double)TMath::ATan2((double)Py,(double)Px);
    if (phi<0)  phi = phi + 6.2831853;
    return phi;
}

// getEta
double CParticle::calcPhiWrapGrad() {
    double  phi=(double)TMath::ATan2((double)Py,(double)Px);
    if (phi<0)  phi = phi + 6.2831853;
    phi = phi* 57.295779578;
    return phi;
}


// getEta
void CParticle::print() {
    printf("CParticle : px=%6.3f py=%6.3f pz=%6.3f cha=%2d m=%6.4f e=%6.3f\n",Px/1000.0,Py/1000.0,Pz/1000.0,Charge,Mass/1000.,E/1000.);

}


// print info about this object
string CParticle::toString(){

      ostringstream outs;
      outs <<   "CParticle  px=" << px()
                 <<  " py=" << py()
                 <<  " pz=" << pz()
                 <<  " e = "<< e()
                 << endl;

     return  outs.str();

} // end of print




// get invariant mass assuming that energy is defined 
double  CParticle::invMassE() {

   double  m2 = ((double)E*(double)E-(double)Px*(double)Px-(double)Py*(double)Py-(double)Pz*(double)Pz);
   double  m=-1;
   if (m2>0) m= TMath::Sqrt(m2);
   return m;

}


// get invariant mass assuming that masses are  defined
double  CParticle::invMassM() {
   double  m2 = ((double)E*(double)E-(double)Px*(double)Px-(double)Py*(double)Py-(double)Pz*(double)Pz);
   double  m=-1;
   if (m2>0) m= TMath::Sqrt(m2);
   return m;

}



