/***************************************************************************
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
#include <TClonesArray.h>

#define MAGENTA "\033[35m"      /* Magenta */
#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */


// custom 
#include "LParticle.h"
#include "CParticle.h"


struct stat sb;


const double kPI   = TMath::Pi();
const double k2PI  = 2*kPI;
using namespace std;

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
         << " expected to provide one input (ROOT) and one output (ROOT) file name and configuration name. \n"
         << " Program stopped! " << endl;
    return 1;
  }

   cout << "ADFilter: Input  File =" << argv[1] << endl;
   cout << "ADFilter: Output File =" << argv[2] << endl;
   cout << "ADFilter: Configuration =" << argv[3] << endl;

	// fastjet
        float CMS=8000;
        const bool   debug=false;
	const float  angNorm=0.15; // factor to increase weights of angles compare to masses
	// determined from analysis of QCD data <energy>/<mass>
	const double ptLepton=30;
	const double ptJet=30;
	const double etaJet=2.4;
	const double etaLep=2.4;
        double PT_LEPTON_LEAD=60; // leading pT of lepton used for selection
        int NN=0;

        int iconfig = atoi( argv[3] );
        //cout << "Configuration name=" <<  iconfig << endl;

	cout << "min PT lepton=" << ptLepton << endl;
	cout << "min PT jet=" << ptJet << endl;
	cout << "min MET   =" << ptJet << endl;
	cout << "max ETA jet=" << etaJet << endl;

        cout << "=================== event selection ================================" << endl;
        if (iconfig == 0)  cout << "  Single lepton with pT>60 GeV as in the ATLAS PRL paper" << endl; 
        if (iconfig == 1)  cout << "  MET>200 GeV trigger" << endl;
        if (iconfig == 2)  cout << "  Single lepton with pT>60 GeV trigger" << endl;
        if (iconfig == 3)  cout << "  2 leptons with pT>30 GeV trigger" << endl;
        if (iconfig == 4)  cout << "  Single photon with pT>150 GeV trigger" << endl;
        if (iconfig == 5)  cout << "  2 photons with pT>30 GeV trigger" << endl;
        if (iconfig == 6)  cout << "  Single jet with pT>500 GeV trigger" << endl;
        if (iconfig == 7)  cout << "  4 jets with pT>200 GeV and 100 GeV for other 3 jets " << endl;
        if (iconfig == 1000)  cout << "  Inclusive jets with pT>50 GeV and 2->2 ME at 1 TeV as in S.Chekanov, R.Zhang, Eur. Phys. J. Plus (2024) 139:237 " << endl;
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

	// get input list 
	std::vector<std::string> files;
        files.push_back( string(argv[1]));

	// output
        TFile * RootFile = new TFile(argv[2], "RECREATE", "Histogram file");
	TH1D * h_debug = new TH1D("EventFlow", "EventFlow", 10, 0, 10.);
        TH1D * h_cross = new TH1D("cross", "Cross [pb], Lumi [pb-1]", 5, 0, 5.);
	TH1D * h_info = new TH1D("Info", "Info", 10, 0, 10.);
	TH2D * h_proj = new TH2D("projection", "projection", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
	TH2D * h_events = new TH2D("events", "events", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
	TProfile2D * h_prof = new TProfile2D("profile", "profile", mSize, 0, (double)(mSize), mSize, 0, (double)mSize, 0, 1000);

	TH1D * h_met = new TH1D("met_pt", "MET",100,0,1000);
        TH1D * h_jetN = new TH1D("jetN", "Nr of light jets",20,0,20);
        TH1D * h_bjetN = new TH1D("bjetN", "Nr of B-Jet",20,0,20);
        TH1D * h_photonN = new TH1D("photonN", "Nr of photons",20,0,20);
        TH1D * h_muonN = new TH1D("muonN", "Nr of muons",20,0,20);
        TH1D * h_electronN = new TH1D("electronN", "Nr of electrons",20,0,20);

        TH1D * h_eta_jet = new TH1D("jet_eta", "jet eta", 40, -5, 5);
        TH1D * h_pt_jet = new TH1D("jet_pt", "pt",100,0,1000);
        TH1D * h_pt_bjet = new TH1D("bjet_pt", "pt",100,0,1000);
        TH1D * h_eta_bjet = new TH1D("bjet_eta", "bjet eta", 40, -5, 5);
        TH1D * h_pt_photon = new TH1D("photon_pt", "leading photon pT",100,0,1000);
        TH1D * h_pt_electron = new TH1D("electron_pt", "leading electron pT",100,0,1000);
        TH1D * h_pt_muon = new TH1D("muon_pt", "leading muon pT",100,0,1000);

         h_pt_jet->GetXaxis()->SetTitle("pT [GeV]");
         h_pt_jet->GetYaxis()->SetTitle("Events");
         h_pt_bjet->GetXaxis()->SetTitle("pT [GeV]");
         h_pt_bjet->GetYaxis()->SetTitle("Events");
	 h_pt_photon->GetXaxis()->SetTitle("pT [GeV]");
         h_pt_photon->GetYaxis()->SetTitle("Events");
         h_pt_electron->GetXaxis()->SetTitle("pT [GeV]");
         h_pt_electron->GetYaxis()->SetTitle("Events");
         h_pt_muon->GetXaxis()->SetTitle("pT [GeV]");
         h_pt_muon->GetYaxis()->SetTitle("Events");
         h_eta_jet->GetXaxis()->SetTitle("Eta");
         h_eta_jet->GetYaxis()->SetTitle("Events");
         

	h_info->Fill("IConfig", iconfig);
        h_info->Fill("PTjet_min", ptJet);
        h_info->Fill("PTlepton_min", ptLepton);
        h_info->Fill("PTgamma_min", ptLepton);

        for (int i=0; i<10; i++)
                  h_info->SetBinError(i,0);


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


    TFile inputFile(argv[1], "READ");
    if (!inputFile.IsOpen()) {
        std::cerr << "Error: Could not open file " <<  argv[1]  << std::endl;
        return 1;
    }


    TH1F* meta = (TH1F*)inputFile.Get("meta");
    CMS = meta->GetBinContent(1); 
    cout << "Read CMS energy [GeV] = " << CMS << endl;


    // Get the TTree from the root file
    TTree* ntuple = dynamic_cast<TTree*>(inputFile.Get("Ntuple"));
    if (!ntuple) {
        std::cerr << "Error: Could not find TTree 'Ntuple' in file " << argv[1]  << std::endl;
        return 1;
    }

    // CM energy in GeV
    Double32_t cmsEnergy;


   // Declaration of leaf types
   Int_t           JET_n;
   vector<Double32_t>  *JET_pt;
   vector<Double32_t>  *JET_eta;
   vector<Double32_t>  *JET_phi;
   vector<Double32_t>  *JET_mass;
   Int_t           bJET_n;
   vector<Double32_t>  *bJET_pt;
   vector<Double32_t>  *bJET_eta;
   vector<Double32_t>  *bJET_phi;
   vector<Double32_t>  *bJET_mass;
   Int_t           EL_n;
   vector<Double32_t>  *EL_pt;
   vector<Double32_t>  *EL_eta;
   vector<Double32_t>  *EL_phi;
   Int_t           MU_n;
   vector<Double32_t>  *MU_pt;
   vector<Double32_t>  *MU_eta;
   vector<Double32_t>  *MU_phi;
   Int_t           PH_n;
   vector<Double32_t>  *PH_pt;
   vector<Double32_t>  *PH_eta;
   vector<Double32_t>  *PH_phi;
   vector<Double32_t>  *PH_e;
   vector<Double32_t>  *MET_eta;
   vector<Double32_t>  *MET_phi;
   vector<Double32_t>  *MET_met;
   vector<Double32_t>  *Evt_Weight;


      // List of branches
   TBranch        *b_JET_n;   //!
   TBranch        *b_JET_pt;   //!
   TBranch        *b_JET_eta;   //!
   TBranch        *b_JET_phi;   //!
   TBranch        *b_JET_mass;   //!
   TBranch        *b_bJET_n;   //!
   TBranch        *b_bJET_pt;   //!
   TBranch        *b_bJET_eta;   //!
   TBranch        *b_bJET_phi;   //!
   TBranch        *b_bJET_mass;   //!
   TBranch        *b_EL_n;   //!
   TBranch        *b_EL_pt;   //!
   TBranch        *b_EL_eta;   //!
   TBranch        *b_EL_phi;   //!
   TBranch        *b_MU_n;   //!
   TBranch        *b_MU_pt;   //!
   TBranch        *b_MU_eta;   //!
   TBranch        *b_MU_phi;   //!
   TBranch        *b_PH_n;   //!
   TBranch        *b_PH_pt;   //!
   TBranch        *b_PH_eta;   //!
   TBranch        *b_PH_phi;   //!
   TBranch        *b_PH_e;   //!
   TBranch        *b_MET_eta;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_met;   //!
   TBranch        *b_Evt_Weight;   //!


   // Set object pointer
   JET_pt = 0;
   JET_eta = 0;
   JET_phi = 0;
   JET_mass = 0;
   bJET_pt = 0;
   bJET_eta = 0;
   bJET_phi = 0;
   bJET_mass = 0;
   EL_pt = 0;
   EL_eta = 0;
   EL_phi = 0;
   MU_pt = 0;
   MU_eta = 0;
   MU_phi = 0;
   PH_pt = 0;
   PH_eta = 0;
   PH_phi = 0;
   PH_e = 0;
   MET_eta = 0;
   MET_phi = 0;
   MET_met = 0;
   Evt_Weight = 0;


   ntuple->SetBranchAddress("JET_n", &JET_n, &b_JET_n);
   ntuple->SetBranchAddress("JET_pt", &JET_pt, &b_JET_pt);
   ntuple->SetBranchAddress("JET_eta", &JET_eta, &b_JET_eta);
   ntuple->SetBranchAddress("JET_phi", &JET_phi, &b_JET_phi);
   ntuple->SetBranchAddress("JET_mass", &JET_mass, &b_JET_mass);
   ntuple->SetBranchAddress("bJET_n", &bJET_n, &b_bJET_n);
   ntuple->SetBranchAddress("bJET_pt", &bJET_pt, &b_bJET_pt);
   ntuple->SetBranchAddress("bJET_eta", &bJET_eta, &b_bJET_eta);
   ntuple->SetBranchAddress("bJET_phi", &bJET_phi, &b_bJET_phi);
   ntuple->SetBranchAddress("bJET_mass", &bJET_mass, &b_bJET_mass);
   ntuple->SetBranchAddress("EL_n", &EL_n, &b_EL_n);
   ntuple->SetBranchAddress("EL_pt", &EL_pt, &b_EL_pt);
   ntuple->SetBranchAddress("EL_eta", &EL_eta, &b_EL_eta);
   ntuple->SetBranchAddress("EL_phi", &EL_phi, &b_EL_phi);
   ntuple->SetBranchAddress("MU_n", &MU_n, &b_MU_n);
   ntuple->SetBranchAddress("MU_pt", &MU_pt, &b_MU_pt);
   ntuple->SetBranchAddress("MU_eta", &MU_eta, &b_MU_eta);
   ntuple->SetBranchAddress("MU_phi", &MU_phi, &b_MU_phi);
   ntuple->SetBranchAddress("PH_n", &PH_n, &b_PH_n);
   ntuple->SetBranchAddress("PH_pt", &PH_pt, &b_PH_pt);
   ntuple->SetBranchAddress("PH_eta", &PH_eta, &b_PH_eta);
   ntuple->SetBranchAddress("PH_phi", &PH_phi, &b_PH_phi);
   ntuple->SetBranchAddress("PH_e", &PH_e, &b_PH_e);
   ntuple->SetBranchAddress("MET_eta", &MET_eta, &b_MET_eta);
   ntuple->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   ntuple->SetBranchAddress("MET_met", &MET_met, &b_MET_met);
   ntuple->SetBranchAddress("Evt_Weight", &Evt_Weight, &b_Evt_Weight);


    bool  err_weight=true;

    float nsum=0;
    float xsum=0;

    // Loop over all entries in the TTree
    Long64_t numEntries = ntuple->GetEntries();
    for (Long64_t i = 0; i < numEntries; ++i) {
        ntuple->GetEntry(i);


        ntot++;
        h_debug->Fill("Input",1);
	
	if (Evt_Weight) {
            for (size_t j = 0; j < Evt_Weight->size(); ++j) {
                //std::cout << "  Weight " << j << ": " << (*weightVector)[j] << endl;
	    }
            if (Evt_Weight->size()>0)  m_weight = Evt_Weight->at(0); 
            
	} else {
                err_weight=true;  
    		//std::cerr << "Error: Weight vector is null." << std::endl;
                m_weight =1.0;
	}

        if (i%1000 == 0)  std::cout << "Event= " << i << " accepted=" <<   NN << " CMS=" << CMS << std::endl;
        vector<LParticle> jets;
        vector<LParticle> bjets; // jets with b-quarks
        vector<LParticle> alljets;
        vector<LParticle> missing;
        vector<LParticle> electrons;
        vector<LParticle> muons;
        vector<LParticle> photons;


	double  met_pt=0;
	double  met_phi=0;
        for (unsigned int j = 0; j < MET_met->size(); ++j) {
                        Double32_t pt=MET_met->at(j);
                        Double32_t eta=0; // MET_eta->at(j);
                        Double32_t phi=MET_phi->at(j);
			// missing
                        TLorentzVector l;
                        met_pt= pt;
                        met_phi = phi;
                        // MET below pT jet is not considered
                        if (met_pt<ptJet){
                                met_pt=0;
                                met_phi=0;
                        }
                        l.SetPtEtaPhiM(met_pt, eta, met_phi, 0);
                        h_met->Fill(met_pt);
                        l.SetPz(0);
                        LParticle p;
                        p.SetP(l);
                        p.SetType(1);
                        missing.push_back(p);
                        break;
	};


        // Print TLorentzVectors for jets
        //std::cout << " Jets:" << std::endl;
	for (int j = 0; j < JET_n; ++j) {
        Double32_t pt=JET_pt->at(j);
	Double32_t eta=JET_eta->at(j);
        Double32_t phi=JET_phi->at(j);
        Double32_t mass=JET_mass->at(j);
	if (iconfig<8  && (pt< ptJet || abs( eta ) > etaJet) ) continue; 
	 // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(pt, eta, phi ,mass);
         h_pt_jet->Fill(pt);
	 h_eta_jet ->Fill(eta); 
	 LParticle p;
         p.SetP(l);
         jets.push_back(p); // light-flavored jets
         alljets.push_back(p);
	 xsum=xsum+(pt);
	 nsum=nsum+1;
	 //cout << pt  << endl;
	}


        for (int j = 0; j < bJET_n; ++j) {
        Double32_t pt=bJET_pt->at(j);
        Double32_t eta=bJET_eta->at(j);
        Double32_t phi=bJET_phi->at(j);
        Double32_t mass=bJET_mass->at(j);
        if (iconfig<8  && (pt< ptJet || abs( eta ) > etaJet) ) continue;
         // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(pt, eta, phi ,mass);
         h_pt_bjet->Fill(pt);
         h_eta_bjet ->Fill(eta);
         LParticle p;
         p.SetP(l);
         bjets.push_back(p); // b jets
         alljets.push_back(p);
         //cout << pt  << endl;
        }




        // Print TLorentzVectors for electrons
        //std::cout << " Electrons:" << std::endl;
        for (int j = 0; j < EL_n; ++j) {
        Double32_t pt=EL_pt->at(j);
        Double32_t eta=EL_eta->at(j);
        Double32_t phi=EL_phi->at(j);
        Double32_t mass=0.000510;
        if (iconfig<8  && (pt< ptLepton || abs( eta  )> etaLep) ) continue; 
        // fill
        TLorentzVector l;
        l.SetPtEtaPhiM(pt,eta, phi, mass);
        h_pt_electron->Fill(pt);
	LParticle p;
        p.SetP(l);
        electrons.push_back(p); // light-flavored jets
	}




       // Print TLorentzVectors for electrons
       //std::cout << " Muons:" << std::endl;
        for (int j = 0; j < MU_n; ++j) {
        Double32_t pt=MU_pt->at(j);
        Double32_t eta=MU_eta->at(j);
        Double32_t phi=MU_phi->at(j);
        Double32_t mass=0.1057;
        if (iconfig<8  && (pt< ptLepton || abs( eta  )> etaLep) ) continue;
        // fill
        TLorentzVector l;
        l.SetPtEtaPhiM(pt,eta, phi, mass);
        h_pt_muon->Fill(pt);
	LParticle p;
        p.SetP(l);
        muons.push_back(p); // light-flavored jets
        }


 
	//std::cout << " Photons:" << std::endl;
        for (int j = 0; j < PH_n; ++j) {
        Double32_t pt=PH_pt->at(j);
        Double32_t eta=PH_eta->at(j);
        Double32_t phi=PH_phi->at(j);
        Double32_t mass=0.0;
        if (iconfig<8  && (pt< ptLepton || abs( eta  )> etaLep) ) continue;
        // fill
        TLorentzVector l;
        l.SetPtEtaPhiM(pt,eta, phi, mass);
        h_pt_photon->Fill(pt);
        LParticle p;
        p.SetP(l);
        photons.push_back(p); // light-flavored jets
        }



	                 // sort all vectors with objects
                        if (jets.size()>1)  std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;
                        if (muons.size()>1) std::sort(muons.begin(), muons.end(), greater<LParticle>() ) ;
                        if (electrons.size()>1) std::sort(electrons.begin(), electrons.end(), greater<LParticle>() ) ;
                        if (photons.size()>1) std::sort(photons.begin(), photons.end(), greater<LParticle>() ) ;
                        if (bjets.size()>1)  std::sort(bjets.begin(), bjets.end(), greater<LParticle>() ) ;
                        if (alljets.size()>1)  std::sort(alljets.begin(), alljets.end(), greater<LParticle>() ) ;

			  h_jetN->Fill(jets.size());
			  h_bjetN->Fill(bjets.size());
			  h_photonN->Fill(photons.size());
			  h_muonN->Fill(muons.size());
			  h_electronN->Fill(electrons.size());

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
                          }

                        if (electrons.size()>0) {
                             LParticle LPP1=electrons.at(0);
                             TLorentzVector L=LPP1.GetP();
                             if (L.Perp()>PT_LEPTON_LEAD) SelectEvent=true;
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
                            if (photons.size() >0) {
                            LParticle p1=photons.at(0);
                            TLorentzVector P1=p1.GetP();
                            if (P1.Perp()>150) SelectEvent=true;
                            };
             	    };



                     // at least 2 photons with pT>30 GeV
                     if (iconfig == 5)
                            if (photons.size() >1) SelectEvent=true;

		     


                 // single jet pT> 500 GeV
                 if (iconfig  == 6) {
                   if (alljets.size() >0) { 
                     LParticle j1=alljets.at(0);
                     TLorentzVector L1=j1.GetP();
                     if (L1.Perp()>500) SelectEvent=true;
                   }
 		   };


                // single jet pT> 30 GeV
                 if (iconfig  == 1000) {
                   if (alljets.size() >0) {
                     LParticle j1=alljets.at(0);
                     TLorentzVector L1=j1.GetP();
                     if (L1.Perp()>30) SelectEvent=true;
                   }
                   };


          // multijets. 1 jet above 200 and other 3 above 100 GeV
       if (iconfig == 7) {
        int II=0;
        if (alljets.size() >3) { 
        LParticle j1=alljets.at(0);
        TLorentzVector L1=j1.GetP();
        if (L1.Perp()>200) II++;

        LParticle j2=alljets.at(1);
        TLorentzVector L2=j2.GetP();
        if (L2.Perp()>100) II++;

        LParticle j3=alljets.at(2);
        TLorentzVector L3=j3.GetP();
        if (L3.Perp()>100) II++;

        LParticle j4=alljets.at(3);
        TLorentzVector L4=j4.GetP();
        if (L4.Perp()>100) II++;
        };

        if (II>3) SelectEvent=true;

       }

        // main cut to select event
        if (SelectEvent == false) continue;



                        h_debug->Fill("Selected",1);

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

                        //cout << met_pt << endl;

			 // we accept non-empty events only!
                        if (met_pt==0 && jets.size()==0 && muons.size()==0 && electrons.size() == 0 && photons.size()==0 && bjets.size()==0) continue;

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


                        if (non_empty>0) h_debug->Fill("Final NonEmpty",1);

                        event++;
                        m_tree->Fill();

        //std::cout << std::endl;
    }

    // Close the input file
    inputFile.Close();

    double xcross=100; // just dummy for now 
    cout << "Average cross section for all files = " << xcross << " pb"<< endl;
    //double width=h_pt->GetBinWidth(1);
    double lumi=(double)ntot/xcross;
    cout << "Lumi for all files  =" << lumi << " pb-1" << endl;
    h_cross->Fill("Cross section [pb]",(double)xcross); // 1
    h_cross->Fill("Lumi [pb-1]",(double)lumi); // 3


    RootFile->Write();
    RootFile->Print();
    RootFile->Close();

     if (err_weight) cout << "Weights for events are not set. Assumes weights=1 for all" << endl;
     cout << "Writing ROOT file " <<  argv[2] << endl;
     cout << "Total events =" << ntot << endl;
     cout << "Selected events =" << event  << endl;


      cout << "Average jet pT=" << xsum/nsum << endl;

      return 0;
}
