/***************************************************************************
 *  How to use ProMC files from HepSim, and how to build anti-KT jets 
 *  S.Chekanov (ANL) chekanov@anl.gov
 *  A library for HEP events storage and processing based on Google's PB   
 *  The project web site: http://atlaswww.hep.anl.gov/hepsim/
****************************************************************************/

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TMath.h"
#include"time.h"
#include <dirent.h>
#include <string>
#include <vector>
#include <map>
#include "TNtupleD.h"
#include <stdio.h>
#include <stdlib.h>
#include <TProfile2D.h>
#include <algorithm>
#include <TRandom2.h>

#define MAGENTA "\033[35m"      /* Magenta */
#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */

struct stat sb;

// promc
#include "ProMC.pb.h"
#include "ProMCBook.h"
#include "LParticle.h"
#include "CParticle.h"

const double kPI   = TMath::Pi();
const double k2PI  = 2*kPI;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

using namespace std;
using namespace promc;

// project event
float**  projectevent(const float  CMS, const int maxN, const int maxNumberTypes,
                      const vector<LParticle> missing,
                      const vector<LParticle> jets,
                      const vector<LParticle> bjets,
                      const vector<LParticle> muons,
                      const vector<LParticle> electrons,
                      const vector<LParticle> photons);


// find all files inside a directory
std::vector<std::string> open(std::string name = "data.in") {
	vector<std::string> ntup;
	ifstream myfile;
	myfile.open(name.c_str(), ios::in);

	if (!myfile) {
		cerr << " -> Can't open input file:  " << name << endl;
		exit(1);
	} else {
		cout << "-> Read data file=" << name << endl;
	}

	string temp;
	while (myfile >> temp) {
		//the following line trims white space from the beginning of the string
		temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), not1(ptr_fun<int, int>(isspace))));
		if (temp.find("#") == 0) continue;
		ntup.push_back(temp);
	}
	cout << "-> Number of files=" << ntup.size()  << endl;
	myfile.close();

	for (unsigned int i=0; i<ntup.size(); i++) {
		cout << ".. file to analyse="+ntup[i] << endl;
	}
	return ntup;
}


// main example
int main(int argc, char **argv)
{



if (argc != 4) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide one input (promc) and one output (ROOT) file name and configuration name. \n"
         << " Program stopped! " << endl;
    return 1;
  }

   cout << "HepSim: Input  File =" << argv[1] << endl;
   cout << "HepSim: Output File =" << argv[2] << endl;
   cout << "HepSim: Configuration =" << argv[3] << endl;

	// fastjet
	const bool   debug=false;
	const float CMS=13000; // GeV;
	//const float CMS=1; // GeV;
	const float  angNorm=0.15; // factor to increase weights of angles compare to masses
	// determined from analysis of QCD data <energy>/<mass>
	const double ptLepton=30;
	const double ptJet=30;
	const double etaJet=2.4;
	const double etaLep=2.4;
	const double Rsize = 0.4;
	const double Rlep_iso=0.1;
	const double Rjet_iso=0.4;
        const double misRateMu=0.1; // 0.1% 
        const double misRateEle=1;  // 1% 
        double PT_LEPTON_LEAD=60; // leading pT of lepton used for selection
        int NN=0;


        int iconfig = atoi( argv[3] );
        //cout << "Configuration name=" <<  iconfig << endl;



       // https://arxiv.org/pdf/1512.01094.pdf
        //1+0.003*pT  (~1000 GeV is about 4%)
        const double misRateB=1;   // 1% b-tag misidentification rate

        TRandom2 misrate;

	cout << "min PT lepton=" << ptLepton << endl;
	cout << "min PT jet=" << ptJet << endl;
	cout << "min MET   =" << ptJet << endl;
	cout << "max ETA jet=" << etaJet << endl;
	cout << "lepton dR isolation=" << Rlep_iso << endl;
	cout << "jet    dR isolation=" << Rjet_iso << endl;
        cout << "Mis-identification rates in %. Ele:" << misRateEle << " Mu:" << misRateMu << endl;
        cout << "mistag b-jet rate % :" << misRateB << endl;

        cout << "=================== event selection ================================" << endl;
        if (iconfig == 0)  cout << "  Single lepton with pT>60 GeV as in the ATLAS PRL paper" << endl; 
        if (iconfig == 1)  cout << "  MET>200 GeV trigger" << endl;
        if (iconfig == 2)  cout << "  Single lepton with pT>60 GeV trigger" << endl;
        if (iconfig == 3)  cout << "  2 leptons with pT>30 GeV trigger" << endl;
        if (iconfig == 4)  cout << "  Single photon with pT>150 GeV trigger" << endl;
        if (iconfig == 5)  cout << "  2 photons with pT>30 GeV trigger" << endl;
        if (iconfig == 6)  cout << "  Single jet with pT>500 GeV trigger" << endl;
        if (iconfig == 7)  cout << "  4 jets with pT>200 GeV and 100 GeV for other 3 jets " << endl;
	cout << "=================== end selection ==================================" << endl;


        int maxNumber=10;   // max number for each object (MET is not counted)
        int maxTypes=5;   // max numbers of types (met not counted)
        // The RMM is 4t3n (maxTypes(t) maxNumber(n)) 
	string names[maxTypes+1] = {"MET","j", "b", "#mu", "e", "#gamma"};
	string names_debug[maxTypes+1] = {"e/met","j", "b", "m", "e", "g"};
	cout << "Project using max number of each object=" << maxNumber << endl;
	cout << "Number of particle types==" << maxTypes << endl;
	const int mSize=maxTypes*maxNumber+1;
	// non-zero in triangular matrix
	int NonZero=((1+mSize)*mSize)/2;


	std::vector<string> Names1;
	Names1.push_back(names[0]);

	for (int h = 1; h < maxTypes+1; h++) {
		for (int i = 1; i <  maxNumber+1; i++) {
			ostringstream ss;
			ss << i;
			Names1.push_back(names[h]+"_{"+ss.str()+"}");
		}
	}


	std::vector<string> Names2;
	for (unsigned int i=0; i<Names1.size(); i++) {
		//cout << "Name=" << i << " " << Names1.at(i) << endl;
		Names2.push_back(Names1.at(i));
	}

	// for plotting
	std::reverse(Names1.begin(), Names1.end());

	int nDiJets=0;
	// total events
	int ntot=0; // total events
	int nfiles=0; // total files
	double weight=1.0;


	//std::vector<std::string> files = open("data.in");

	// get input list 
	std::vector<std::string> files;
        files.push_back( string(argv[1]));

	double cross=0;
	double xcross=0;
	int nselect=0;

	//string outputfile="output.root";
	//cout << "\n -> Output file is =" << outputfile << endl;
	//TFile * RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
	//  RootFile->SetCompressionLevel(0);
        TFile * RootFile = new TFile(argv[2], "RECREATE", "Histogram file");

	TH1D * h_debug = new TH1D("EventFlow", "EventFlow", 10, 0, 10.);
	TH1D * h_cross = new TH1D("cross", "cross,events,lumi", 5, 0, 5.);
	TH1D * h_info = new TH1D("dijet_info", "Dijet info", 5, 0, 5.);
	TH1D * h_n_lepton = new TH1D("lepton_nr", "Nr of leptons",20,0,20);
	TH1D * h_iso = new TH1D("iso_energy", "isolation fraction", 100, 0, 3.0);
	TH2D * h_proj = new TH2D("projection", "projection", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
	TH2D * h_events = new TH2D("events", "events", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
	TProfile2D * h_prof = new TProfile2D("profile", "profile", mSize, 0, (double)(mSize), mSize, 0, (double)mSize, 0, 1000);

	TH1D * h_met = new TH1D("met_pt", "MET",100,0,1000);
        TH1D * h_jetN = new TH1D("jetN", "Nr of light jets",20,0,20);
        TH1D * h_bjetN = new TH1D("bjetN", "Nr of B-Jet",20,0,20);
        TH1D * h_photonN = new TH1D("photonN", "Nr of photons",20,0,20);
        TH1D * h_muonN = new TH1D("muonN", "Nr of muons",20,0,20);
        TH1D * h_electronN = new TH1D("electronN", "Nr of electrons",20,0,20);

        TH1D * h_eta_jet = new TH1D("jet_eta", "eta", 40, -10, 10);
        TH1D * h_pt_jet = new TH1D("jet_pt", "pt",100,0,1000);
        TH1D * h_pt_bjet = new TH1D("bjet_pt", "pt",100,0,1000);
        TH1D * h_eta_bjet = new TH1D("bjet_eta", "eta", 40, -10, 10);
        TH1D * h_pt_photon = new TH1D("photon_pt", "leading photon pT",100,0,1000);
        TH1D * h_pt_electron = new TH1D("electron_pt", "leading electron pT",100,0,1000);
        TH1D * h_pt_muon = new TH1D("muon_pt", "leading muon pT",100,0,1000);



        // remember which matrix do you fill
        TH1D * h_dimensions = new TH1D("dimensions", "(1)maxNumber,(2)maxTypes, (3)mSize",5,0,5);
        h_dimensions->Fill(1,(float)maxNumber);
        h_dimensions->Fill(2,(float)maxTypes);
        h_dimensions->Fill(3,(float)mSize);

	// each accepted event
	int event=0;
	const int nevhisto=50;
	TH2D * h_proje[ nevhisto];
	for(int i=0; i<nevhisto; i++)
		h_proje[i] = new TH2D(Form("event_%02d",i),Form("event_%02d",i), mSize, 0, (double)(mSize), mSize, 0, (double)mSize);



	// initialize with 0
	for (int h = 0; h < mSize; h++) {
		for (int w = 0; w < mSize; w++) {
			int i1=h;
			int i2=w;
			h_proj->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),0);
			h_prof->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),0);
			h_events->Fill((Names2.at(i1)).c_str(), (Names1.at(i2)).c_str(),0);
			for(int i=0; i<nevhisto; i++) h_proje[i]->Fill((Names2.at(i1)).c_str(), (Names1.at(i2)).c_str(),0);
		}}




	TTree *  m_tree  = new TTree("inputNN","inputNN");
	m_tree->SetAutoFlush(100000);
	Int_t m_id;
	std::vector<Double32_t> m_proj;
        std::vector<UInt_t> m_proj_index1;
        std::vector<UInt_t> m_proj_index2;
        Int_t m_run;
        Int_t m_event;
        Double32_t m_weight;

	std::vector<Double32_t> m_multi;
	m_tree->Branch("id",  &m_id);
        m_tree->Branch("run",  &m_run);
        m_tree->Branch("event",  &m_event);
        m_tree->Branch("weight",  &m_weight);
        m_tree->Branch("proj",   &m_proj);
        m_tree->Branch("proj_index1",   &m_proj_index1);
        m_tree->Branch("proj_index2",   &m_proj_index2);
         
	m_tree->Branch("multiplicity",   &m_multi);

	// jets
	Strategy strategy = fastjet::Best;
	JetDefinition jet_def(fastjet::antikt_algorithm, Rsize, strategy);
	//JetDefinition jet_def(fastjet::kt_algorithm, Rparam, strategy);
	//JetDefinition jet_def(fastjet::cambridge_algorithm, Rparam, strategy);
	//JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);


	std::map<int,int> chargemap; // pad that keeps charge*3

	for(unsigned int m=0; m < files.size(); m++){
		string Rfile=files[m];
		ProMCBook*  epbook = new ProMCBook(Rfile.c_str(),"r");

		cout << "\n\n Start to read.." << endl;
		// get the version number
		int  h=epbook->getVersion();
		if (m==0) cout << "Version = " << h << endl;
		// get the description of this file
		string d=epbook->getDescription();
		if (m==0) cout << "Description = " << d << endl;
		int  nev=epbook->getEvents();
		cout << "Events = " << nev  << endl;
		// get the header file with units, cross sections etc.
		ProMCHeader header = epbook->getHeader();


		// this part reads the header information with particle data.
		// you can access names, id, charge*3, and other useful
		// information from PDG table. As an example, we make a map
		// that keeps charge*3.
		if (m==0) {
			for (int i1=0; i1<header.particledata_size(); i1++){
				ProMCHeader_ParticleData data= header.particledata(i1);
				int charge3=data.charge();
				int id=data.id();
				double mass=data.mass();
				string name=data.name();
				//cout << "name=" << name << " mass=" << mass << " charge=" << charge3 << endl;
				chargemap[id]=charge3;
			}
		}


		// here are the units
		double kEV=(double)(header.momentumunit());
		//double kLe=(double)(header.lengthunit());


		// loop over all events
		for (int nn=0; nn<nev; nn++){
			if (epbook->next() !=0) continue;
			ProMCEvent eve = epbook->get();


			// get truth information
			ProMCEvent_Particles  *pa=eve.mutable_particles();
			h_debug->Fill("Events",1.0);
			ntot++;

			double xlumi=0;
			if (m>0) {
				xlumi=(double)ntot/cross; // lumi so far
			}

			if (ntot%1000==0)
				cout <<  " # Events=" << ntot  << " X-cross="<< cross << " pb" << " Lumi so far=" << xlumi/1000.0 << " fb-1" << " selected="<< NN << endl;


			vector<PseudoJet> avec;
			vector<int> nrid;
			vector<LParticle> candidates;
                        vector<LParticle> bquarks;

			double pxsum=0;
			double pysum=0;
			double pzsum=0;

			// fill stable and no neutrino
			for (int i2=0; i2<pa->pdg_id_size(); i2++){


				int type=pa->pdg_id(i2);
				if (abs(type)==12 || abs(type)==14 || abs(type)==16 ) continue;
				double px= pa->px(i2)/kEV;
				double py= pa->py(i2)/kEV;
				double pz= pa->pz(i2)/kEV;
				double ee= pa->energy(i2)/kEV;
				double pt=sqrt(px*px+py*py);
				double eta=-log(tan(atan2(pt,(double)pz)/2));


                        // b-quarks
                        if (abs(type) ==5 and pt>0.3*ptJet) {
                                int charge=1;
                                if (type<0) charge=-1;
                                LParticle b(px,py,pz,ee,charge);
                                b.SetStatus(0);
                                b.SetType(type);
                                bquarks.push_back(b);
                        };

                                if (pa->status(i2)!=1) continue;
                                // muon or electron or photon
				if (abs(type) ==11 || abs(type) ==13 || abs(type) ==22 ) {
					int charge=1;
					if (type<0) charge=-1;
					LParticle p(px,py,pz,ee,charge);
					p.SetStatus(i2);
					p.SetType(type);
					if (pt>ptLepton && TMath::Abs(eta)<etaLep) {
						if (debug) cout << i2 << " ProMC type=" << type << " ProMC pt=" << pt << "ProMC eta=" << eta << endl;
						candidates.push_back(p);
					}

				};
				// int charge=chargemap[pa->pdg_id(j)]; // get charge
				if ( pt < 0.1)                   continue;
				if ( fabs(eta)> 3.5 )            continue;
				avec.push_back(  PseudoJet(px,py,pz,ee)  );
				nrid.push_back(i2);

				pxsum=pxsum+px;
				pysum=pysum+py;
				pzsum=pzsum+pz;

			}

			// isolate candidates
			vector<LParticle> leptons;
			double IsoEnergy=0.1;
			// isolate leading lepton
			for (unsigned int k1 = 0; k1<candidates.size(); k1++) {
				LParticle xlep=(LParticle)candidates.at(k1);
				TLorentzVector lep=xlep.GetP();
				double p_pt=lep.Et();
				double p_eta=lep.PseudoRapidity();
				double p_phi=lep.Phi();
				if (p_phi<0) p_phi=k2PI+p_phi;
				if (p_pt<ptLepton) continue;
				double esumP=0;
				for (unsigned int k2 = 0; k2<avec.size(); k2++) {
					PseudoJet part = avec.at(k2);
					double pt=part.perp();
					double eta=part.pseudorapidity();
					double phi=part.phi();
					if (phi<0) phi=k2PI+phi;
					double deta    = p_eta - eta;
					double dphi    = p_phi - phi;
					double adphi=TMath::Abs(dphi);
					if (adphi>kPI) adphi=k2PI-adphi;
					double ddr = TMath::Sqrt(deta*deta+adphi*adphi);
					if (ddr<Rlep_iso) esumP=esumP+pt;
				}
				double isoFrac=esumP/p_pt;
				//cout << "Esum=" << esumP << " p_pt = " << p_pt << " iso=" << isoFrac << endl;
				h_iso->Fill(isoFrac);
				if (isoFrac<1.0+IsoEnergy) {
					leptons.push_back(xlep);
				};
			}


			unsigned int nLeptons=leptons.size();
			if (nLeptons>1) std::sort(leptons.begin(), leptons.end(), greater<LParticle>() ) ;
			h_n_lepton->Fill(nLeptons);


			// remove isolated lepton from vectors used for jets
			vector<PseudoJet> hadrons;
			for (unsigned int k = 0; k<avec.size(); k++) {
				PseudoJet part = avec.at(k);
				int id=nrid.at(k);
				int isLep=false;
				for (unsigned int ll=0; ll<leptons.size(); ll++){
					LParticle LL=leptons.at(ll);
					int id_lep=LL.GetStatus();
					if (id_lep == id) isLep=true;
				}
				if (!isLep) hadrons.push_back(part);
			}




			// make jets
			ClusterSequence clust_seq(hadrons, jet_def);
			vector<PseudoJet> jets_truth = clust_seq.inclusive_jets(ptJet);
			vector<PseudoJet> sorted_jets = sorted_by_pt(jets_truth);
			vector<LParticle> jets;
                        vector<LParticle> bjets; // jets with b-quarks
                        vector<LParticle> alljets;

			for (unsigned int k = 0; k<sorted_jets.size(); k++) {
				double eta=sorted_jets[k].pseudorapidity();
				if ( fabs(eta)> etaJet )            continue;
				double phi=sorted_jets[k].phi();
				double pt = sorted_jets[k].perp();
				if (pt<ptJet)                  continue;
				double e = sorted_jets[k].e();


                               // misidentified b-jets using misRateB rate
                                int FakeB=0;
                                if (100.0*misrate.Rndm()< (misRateB+0.003*pt) ) FakeB=1;

                                // find and label b-quark jets
                                int matchB=0;
                                for (unsigned int k3=0; k3<bquarks.size(); k3++) {
                                        LParticle p=(LParticle)bquarks.at(k3);
                                        TLorentzVector L= p.GetP();
                                        double pt_b = L.Perp();
                                        double eta_b = L.PseudoRapidity();
                                        double phi_b = L.Phi();
                                        if (phi_b<0) phi_b=k2PI+phi_b;
                                        double deta=eta_b-eta;
                                        double dphi=phi_b-phi;
                                        double dR=sqrt(deta*deta + dphi*dphi);
                                        if (dR<Rsize && pt_b>0.5*pt) matchB=1;
                                }


                               // light flavor jets
                                TLorentzVector l;
                                l.SetPtEtaPhiE(pt, eta, phi, e);
                                LParticle p;
                                p.SetP(l);
                                p.SetType(matchB); // label b-jet quarks
                                int nmulti=sorted_jets[k].constituents().size();
                                p.SetCharge(nmulti); // multiplicity
                                alljets.push_back(p);  // all jets

				if (matchB==0 && FakeB==0) {
                                        h_pt_jet->Fill(pt);
                                        h_eta_jet->Fill(eta);
                                        jets.push_back(p); // light-flavored jets
                                } else if (matchB==1 || FakeB==1) {
                                        h_pt_bjet->Fill(pt);
                                        h_eta_bjet->Fill(eta);
                                        if (FakeB==1) p.SetType(-1);
                                        bjets.push_back(p); // b-jets and fake b-jets 
                                }


			} // end loop


			unsigned int nJets=jets.size();
			if (nJets>1) std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;


			// overlap removal for 3 first jets (how many tipes)
			int nMaxJet=jets.size();
			if (nMaxJet>maxNumber) nMaxJet=maxNumber;

			vector<LParticle> leptons_iso;
			for (unsigned int ll=0; ll<nLeptons; ll++){
				LParticle LL=leptons.at(ll);
				TLorentzVector LP=LL.GetP();
				double phi_lep=LP.Phi();
				double eta_lep=LP.PseudoRapidity();
				double y_lep=LP.Rapidity();

				bool found=false; // isolated from 3 leading jets
				for (int ii=0; ii<nMaxJet; ii++){
					LParticle LPP=jets.at(ii);
					TLorentzVector LP=LPP.GetP();
					double phi_jet=LP.Phi();
					double eta_jet=LP.PseudoRapidity();
					double y_jet=LP.Rapidity();
					double deta=TMath::Abs(y_jet-y_lep);
					double dphi=TMath::Abs(phi_jet - phi_lep);
					if (dphi>kPI) dphi=k2PI-dphi;
					double dR=TMath::Sqrt(deta*deta+dphi*dphi);
					if (dR<Rjet_iso) found=true;
				}
				if (found) continue;
				leptons_iso.push_back(LL);
			}


			// now classify isolated candidates
			vector<LParticle> missing;
			vector<LParticle> electrons;
			vector<LParticle> muons;
			vector<LParticle> photons;
			for (unsigned int k1 = 0; k1<leptons_iso.size(); k1++) {
				LParticle xlep=(LParticle)leptons_iso.at(k1);
				int type=xlep.GetType();
				if (abs(type) ==11) electrons.push_back(xlep);
				if (abs(type) ==13) muons.push_back(xlep);
				if (abs(type) ==22) photons.push_back(xlep);
			}



			// missing
                        TLorentzVector l;
                        double met_pt=sqrt(pxsum*pxsum+pysum*pysum);
                        double met_phi = 0;
                        if (pxsum>0)
                                met_phi=atan2(pysum,pxsum);
                        if (pxsum<0)
                                met_phi=atan2(pysum,pxsum) + kPI;
                        // MET below pT jet is not considered
                        if (met_pt<ptJet){
                                met_pt=0;
                                met_phi=0;
                        }

                        l.SetPtEtaPhiM(met_pt, 0, met_phi, 0);
                        h_met->Fill(met_pt);
                        l.SetPz(0);
                        LParticle p;
                        p.SetP(l);
                        p.SetType(1);
                        missing.push_back(p);



                       // sort all vectors with objects 
                        if (jets.size()>1)  std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;
                        if (muons.size()>1) std::sort(muons.begin(), muons.end(), greater<LParticle>() ) ;
                        if (electrons.size()>1) std::sort(electrons.begin(), electrons.end(), greater<LParticle>() ) ;
                        if (photons.size()>1) std::sort(photons.begin(), photons.end(), greater<LParticle>() ) ;
                        if (bjets.size()>1)  std::sort(bjets.begin(), bjets.end(), greater<LParticle>() ) ;
                        if (alljets.size()>1)  std::sort(alljets.begin(), alljets.end(), greater<LParticle>() ) ;



		        /************************************************
			 * Event selection
			 ************************************************/

                        bool SelectEvent=false;


                        // single lepton trigger pT>60 GeV
                        // PRL or Trigger 2
			if (iconfig == 0 || iconfig == 2) {
                        if (muons.size()>0) {
                             LParticle LPP1=muons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             if (L.Perp()>PT_LEPTON_LEAD) SelectEvent=true; 
                             h_pt_muon->Fill(L.Perp());
                          }

                        if (electrons.size()>0) {
                             LParticle LPP1=electrons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             if (L.Perp()>PT_LEPTON_LEAD) SelectEvent=true; 
                             h_pt_electron->Fill(L.Perp());
                          }
                        } // iconfig == 0 



		       // MET trigger
                       if (iconfig == 1) {
                               if (met_pt>200) SelectEvent=true; 
		       } // iconfig == 0 



		        //  2 lepton trigger
                       if (iconfig == 3) {
                       PT_LEPTON_LEAD=25.0;
                       int NrLepton=0;
                       for (unsigned int j1=0; j1<muons.size(); j1++) {
                          LParticle muon=muons.at(j1);
                          TLorentzVector LM=muon.GetP();
                          if (LM.Perp()>PT_LEPTON_LEAD) NrLepton++;
                        }
                        for (unsigned int j1=0; j1<electrons.size(); j1++) {
                           LParticle ele=electrons.at(j1);
                           TLorentzVector LM=ele.GetP();
                           if (LM.Perp()>PT_LEPTON_LEAD)  NrLepton++;
                        }
                           if (NrLepton>1) SelectEvent=true;   
                       };


		      // single photon > 150 GeV
                     if (iconfig == 4) {
                            if (photons.size() < 1) continue;
                            LParticle p1=photons.at(0);
                            TLorentzVector P1=p1.GetP();
                           if (P1.Perp()>150) SelectEvent=true; 
                        };



	             // at least 2 photons with pT>30 GeV
                     if (iconfig == 5)
                            if (photons.size() >1) SelectEvent=true;


	         // single jet pT> 500 GeV 
                 if (iconfig  == 6) {
                    if (alljets.size() < 1) continue;
                   LParticle j1=alljets.at(0);
                   TLorentzVector L1=j1.GetP();
                   if (L1.Perp()<500) continue;
                 };



	  // multijets. 1 jet above 200 and other 3 above 100 GeV
       if (iconfig == 7) {
        if (alljets.size() < 4) continue;
        LParticle j1=alljets.at(0);
        TLorentzVector L1=j1.GetP();
        if (L1.Perp()<200) continue;

        LParticle j2=alljets.at(1);
        TLorentzVector L2=j2.GetP();
        if (L2.Perp()<100) continue;

        LParticle j3=alljets.at(2);
        TLorentzVector L3=j3.GetP();
        if (L3.Perp()<100) continue;

        LParticle j4=alljets.at(3);
        TLorentzVector L4=j4.GetP();
        if (L4.Perp()<100) continue;

    }

             // main cut to select event 
             if (SelectEvent == false) continue;



                        if (photons.size()>0) {
                             LParticle LPP1=photons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             h_pt_photon->Fill(L.Perp());
                          }


                         h_debug->Fill("Final",1.0);

                        NN++;
                        
                        h_jetN->Fill((float)jets.size());
                        h_bjetN->Fill((float)bjets.size());
                        h_photonN->Fill((float)photons.size());
                        h_muonN->Fill((float)muons.size());
                        h_electronN->Fill((float)electrons.size());



			m_proj.clear();
                        m_proj_index1.clear();
                        m_proj_index2.clear();

			m_multi.clear();
			//m_m=muons.size();
			//m_e=electrons.size();
			//m_g=photons.size();
			//m_j=jets.size();
			m_id=0;
                        m_event=NN;
                        m_run=0;
                        m_weight=1.0;



			// we accept non-empty events only!
			if (met_pt==0 && jets.size()==0 && muons.size()==0 && electrons.size() == 0 && photons.size()==0) continue;

			if (event<nevhisto){
				cout << event << " Nr of gamma= " << photons.size() << " e=" << electrons.size() << " mu=" << muons.size() << " jet=" << jets.size() << endl;
			}

			// matrix has size:
			if (debug) {

				cout << "# Nr jets=" << jets.size()<< " muons=" <<  muons.size() << " ele=" << electrons.size() << " pho=" << photons.size() << endl;

				cout << " Nr of Photons=" << photons.size() << endl;
				for (unsigned int i1=0; i1<photons.size(); i1++){
					LParticle LPP1=photons.at(i1);
					TLorentzVector LP1=LPP1.GetP();
					cout << i1<< " pt=" << LP1.Et() << " eta=" << LP1.Eta() << " energy=" << LP1.E() << endl;
				}

				cout << "\n Nr of Jets=" << jets.size() << endl;
				for (unsigned int i1=0; i1<jets.size(); i1++){
					LParticle LPP1=jets.at(i1);
					TLorentzVector LP1=LPP1.GetP();
					cout << i1 << " pt=" << LP1.Et() << " eta=" << LP1.Eta() << " energy=" << LP1.E() << endl;
				}




			};


                        // return rapidity-mass matrix, RMM,  (with zeros)
                        float** projArray =  projectevent(CMS, maxNumber, maxTypes, missing, jets, bjets, muons, electrons, photons);


			// triangular matrix is a special kind of square matrix.
			// upper triangular matrix or right triangular matrix.
			// non-zero in triangular matrix
			if (debug) {
				cout << "  " << names_debug[0] << "          ";
				for (int h = 1; h < maxTypes+1; h++)
					for (int i = 0; i <  maxNumber; i++)
						cout << names_debug[h] << i << "       ";
				cout << endl;
				for (int h = 0; h < mSize; h++) {
					cout << h << " ";
					for (int w = 0; w < mSize; w++)
					{
						if (h != w) printf("%.2e ", float(projArray[h][w]));
						else  {
							//cout << std::setprecision(1) << std::scientific <<  float(projArray[h][w]);
							//printf("%.2e ", float(projArray[h][w]));
							const char* tt="ok";
							cout << RED; printf("%.2e ", float(projArray[h][w])); cout << RESET;
							// else printf("\033[1;31%.2e033[0m", float(projArray[h][w]));
						}

					}
					printf("\n");
				}
			}; // end debug


			int nMax=mSize *mSize;
			Double32_t* array = new Double32_t[nMax]; // allocate
			for (int h = 0; h < nMax; h++) array[h]=0;


			Double32_t* arrayM = new Double32_t[5]; // allocate
			for (int h = 0; h < 5; h++) arrayM[h]=0;

			// simple count of multiplicities:
			// met, Nr(jets), Nr(mu), Nr(el), Nr(gam)
			Double32_t* arrayMuliplicity = new Double32_t[5]; // allocate
			arrayM[0]=met_pt/CMS;
			unsigned int nj=maxNumber;
			if (jets.size()<nj) nj=jets.size();
			arrayM[1]=(Double32_t)(nj/(float)maxNumber);
			nj=maxNumber;
			if (muons.size()<nj) nj=muons.size();
			arrayM[2]=(Double32_t)(nj/(float)maxNumber);
			nj=maxNumber;
			if (electrons.size()<nj) nj=electrons.size();
			arrayM[3]=(Double32_t)(nj/(float)maxNumber);
			nj=maxNumber;
			if (photons.size()<nj) nj=electrons.size();
			arrayM[4]=(Double32_t)(nj/(float)maxNumber);
			for (int i = 0; i < 5; i++) m_multi.push_back( arrayM[i] );


			int count=0;
                        int non_empty=0;

			// only non-zero events 
			for (int h = 0; h < mSize; h++) {
				for (int w = 0; w < mSize; w++) {
					float dd=projArray[w][h];

					if (dd != 0) {
						//if (h<w) h_matAngle->Fill(dd);
						//if (h>w) h_matMasses->Fill(dd);
						//if (h==w) h_matPt->Fill(dd);
						// reduce angle weights
						if (h<w)  {
							dd=dd*angNorm; // decrease angles by 0.01
							//h_matAngNorm->Fill(dd);
						};
					}

					int i1=h;
					int i2=mSize-w-1;

					// profile grabs 0, so average for masses is smaller
					h_prof->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),dd);

					if (event<nevhisto)
						h_proje[event]->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),dd);


					// to get correct masses and suppress 0
					if (dd != 0) {
						h_proj->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),dd);
						h_events->Fill((Names2.at(i1)).c_str(), (Names1.at(i2)).c_str(),1.0);

                                                // prepare for output. Fill non-zero values, row and column index 
                                                m_proj.push_back(dd);
                                                m_proj_index1.push_back( w );
                                                m_proj_index2.push_back( h );
                                                non_empty++;

					}

					array[count]=dd;
					count++;
					//cout << i << " " <<  array[i] << endl;
				}
			}

                        //fill this matrix if you want to store zeros too!
			//for (int i = 0; i < nMax; i++) m_proj.push_back( array[i] );

			delete [] arrayM;
			delete [] array;
                        for (int w1 = 0; w1 < mSize; w1++)  delete[] projArray[w1];
                        delete[] projArray;


                        if (non_empty>0) h_debug->Fill("Final NonEmpty",1.0);

			event++;
			m_tree->Fill();

		} // end event loop


		ProMCStat stat = epbook->getStatistics();
		cross=stat.cross_section_accumulated();
		epbook->close(); // close
		nfiles++;
		xcross=xcross+cross;


	} // end loop over all files



	xcross=xcross/(double)nfiles; // average cross for all files
	cout << "Total files=" << nfiles << endl;
	cout << "Total events =" << ntot << endl;
        cout << "Selected RMM events =" << event  << endl;
	cout << "Average cross section for all files = " << xcross << " pb"<< endl;
	//double width=h_pt->GetBinWidth(1);
	double lumi=(double)ntot/xcross;
	cout << "Lumi for all files  =" << lumi << " pb-1" << endl;
	h_cross->Fill("Cross section",(double)xcross); // 1
	h_cross->Fill("Input Events",(double)ntot); // 2
        h_cross->Fill("Selected Events",(double)event); // 2
	h_cross->Fill("Lumi [pb]",(double)lumi); // 3
	h_cross->Fill("Input Files",(double)nfiles);


	RootFile->Write();
	RootFile->Print();
	RootFile->Close();

	cout << "Writing ROOT file " <<  argv[2] << endl;

	return 0;
}
