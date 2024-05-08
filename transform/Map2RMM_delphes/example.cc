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
	TH1D * h_cross = new TH1D("cross", "cross,events,lumi", 5, 0, 5.);
	TH1D * h_info = new TH1D("dijet_info", "Dijet info", 5, 0, 5.);
	TH1D * h_n_lepton = new TH1D("lepton_nr", "Nr of leptons",20,0,20);
	TH2D * h_proj = new TH2D("projection", "projection", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
	TH2D * h_events = new TH2D("events", "events", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
	TProfile2D * h_prof = new TProfile2D("profile", "profile", mSize, 0, (double)(mSize), mSize, 0, (double)mSize, 0, 1000);

	TH1D * h_met = new TH1D("met_pt", "MET",100,0,1000);
        TH1D * h_jetN = new TH1D("jetN", "Nr of light jets",20,0,20);
        TH1D * h_bjetN = new TH1D("bjetN", "Nr of B-Jet",20,0,20);
        TH1D * h_photonN = new TH1D("photonN", "Nr of photons",20,0,20);
        TH1D * h_muonN = new TH1D("muonN", "Nr of muons",20,0,20);
        TH1D * h_electronN = new TH1D("electronN", "Nr of electrons",20,0,20);

        TH1D * h_eta_jet = new TH1D("jet_eta", "jet eta", 40, -10, 10);
        TH1D * h_pt_jet = new TH1D("jet_pt", "pt",100,0,1000);
        TH1D * h_pt_bjet = new TH1D("bjet_pt", "pt",100,0,1000);
        TH1D * h_eta_bjet = new TH1D("bjet_eta", "bjet eta", 40, -10, 10);
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


    TFile inputFile(argv[1], "READ");
    if (!inputFile.IsOpen()) {
        std::cerr << "Error: Could not open file " <<  argv[1]  << std::endl;
        return 1;
    }


    // Get the TTree from the root file
    TTree* ntuple = dynamic_cast<TTree*>(inputFile.Get("Ntuple"));
    if (!ntuple) {
        std::cerr << "Error: Could not find TTree 'Ntuple' in file " << argv[1]  << std::endl;
        return 1;
    }

    // CM energy in GeV
    Double32_t cmsEnergy;
    // Create TClonesArrays to store TLorentzVectors
    TClonesArray* jetArray = nullptr;
    TClonesArray* bjetArray = nullptr;
    TClonesArray* electronArray = nullptr;
    TClonesArray* muonArray = nullptr;
    TClonesArray* photonArray = nullptr;
    TClonesArray* metArray = nullptr;
    std::vector<Double_t>* weightVector = nullptr;

    // Set branch addresses for each TClonesArray
    ntuple->SetBranchAddress("Jets", &jetArray);
    ntuple->SetBranchAddress("BJets", &bjetArray);
    ntuple->SetBranchAddress("Electrons", &electronArray);
    ntuple->SetBranchAddress("Muons", &muonArray);
    ntuple->SetBranchAddress("Photons", &photonArray);
    ntuple->SetBranchAddress("MET", &metArray);
    ntuple->SetBranchAddress("Weight", &weightVector);
    ntuple->SetBranchAddress("CMS_Energy", &cmsEnergy);


    bool  err_weight=true;

    // Loop over all entries in the TTree
    Long64_t numEntries = ntuple->GetEntries();
    for (Long64_t i = 0; i < numEntries; ++i) {
        ntuple->GetEntry(i);


        ntot++;
	
	if (weightVector) {
            for (size_t j = 0; j < weightVector->size(); ++j) {
                //std::cout << "  Weight " << j << ": " << (*weightVector)[j] << endl;
	    }
            if (weightVector->size()>0)  m_weight = (*weightVector)[0]; 
            
	} else {
                err_weight=true;  
    		//std::cerr << "Error: Weight vector is null." << std::endl;
                m_weight =1.0;
	}


	CMS=(float)cmsEnergy;
        if (i%1000 == 0)  std::cout << "Event= " << i << " accepted=" <<   NN << " CMS=" << CMS << std::endl;
        vector<LParticle> jets;
        vector<LParticle> bjets; // jets with b-quarks
        vector<LParticle> alljets;
        vector<LParticle> missing;
        vector<LParticle> electrons;
        vector<LParticle> muons;
        vector<LParticle> photons;

	               
	                TLorentzVector* met = dynamic_cast<TLorentzVector*>(metArray->At(0));
                        // missing
                        TLorentzVector l;
                        double met_pt= met->Pt(); 
                        double met_phi = met->Phi();
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

        // Print TLorentzVectors for jets
        //std::cout << " Jets:" << std::endl;
        for (int j = 0; j < jetArray->GetEntries(); ++j) {
            TLorentzVector* jet = dynamic_cast<TLorentzVector*>(jetArray->At(j));
            //std::cout << "   Jet " << j << ": (Pt, Eta, Phi, M) = (" << jet->Pt() << ", " << jet->Eta() << ", " << jet->Phi() << ", " << jet->M() << ")" << std::endl;

	 // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(jet->Pt(), jet->Eta(), jet->Phi() , jet->M());
         h_pt_jet->Fill(jet->Pt());
	 h_eta_jet ->Fill(jet->Eta()); 
	 LParticle p;
         p.SetP(l);
         jets.push_back(p); // light-flavored jets
         alljets.push_back(p);

	}

        // Print TLorentzVectors for b-jets
        //std::cout << "  B-Jets:" << std::endl;
        for (int j = 0; j < bjetArray->GetEntries(); ++j) {
            TLorentzVector* bjet = dynamic_cast<TLorentzVector*>(bjetArray->At(j));
            //std::cout << "   B-Jet " << j << ": (Pt, Eta, Phi, M) = (" << bjet->Pt() << ", " << bjet->Eta() << ", " << bjet->Phi() << ", " << bjet->M() << ")" << std::endl;

         // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(bjet->Pt(), bjet->Eta(), bjet->Phi() , bjet->M());
         h_pt_bjet->Fill(bjet->Pt());
	 h_eta_bjet ->Fill(bjet->Eta()); 
	 LParticle p;
         p.SetP(l);
         bjets.push_back(p); // b jets
         alljets.push_back(p);

	}

        // Print TLorentzVectors for electrons
        //std::cout << " Electrons:" << std::endl;
        for (int j = 0; j < electronArray->GetEntries(); ++j) {
            TLorentzVector* electron = dynamic_cast<TLorentzVector*>(electronArray->At(j));
            //std::cout << "   Electron " << j << ": (Pt, Eta, Phi, M) = (" << electron->Pt() << ", " << electron->Eta() << ", " << electron->Phi() << ", " << electron->M() << ")" << std::endl;

         // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(electron->Pt(), electron->Eta(), electron->Phi() , 0.000510);
         h_pt_electron->Fill(electron->Pt());
	 LParticle p;
         p.SetP(l);
         electrons.push_back(p); // light-flavored jets

	}

        // Print TLorentzVectors for muons
        //std::cout << " Muons:" << std::endl;
        for (int j = 0; j < muonArray->GetEntries(); ++j) {
            TLorentzVector* muon = dynamic_cast<TLorentzVector*>(muonArray->At(j));
            //std::cout << "   Muon " << j << ": (Pt, Eta, Phi, M) = (" << muon->Pt() << ", " << muon->Eta() << ", " << muon->Phi() << ", " << muon->M() << ")" << std::endl;
         // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(muon->Pt(), muon->Eta(), muon->Phi(),0.1057);
         h_pt_muon->Fill(muon->Pt());
	 LParticle p;
         p.SetP(l);
         muons.push_back(p); // light-flavored jets

	}

        // Print TLorentzVectors for photons
        //std::cout << " Photons:" << std::endl;
        for (int j = 0; j < photonArray->GetEntries(); ++j) {
            TLorentzVector* photon = dynamic_cast<TLorentzVector*>(photonArray->At(j));
            //std::cout << "   Photon " << j << ": (Pt, Eta, Phi, M) = (" << photon->Pt() << ", " << photon->Eta() << ", " << photon->Phi() << ", " << photon->M() << ")" << std::endl;
         // fill
         TLorentzVector l;
         l.SetPtEtaPhiM(photon->Pt(), photon->Eta(), photon->Phi(),0.0);
         h_pt_photon->Fill(photon->Pt());
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

        //std::cout << std::endl;
    }

    // Close the input file
    inputFile.Close();


     RootFile->Write();
     RootFile->Print();
     RootFile->Close();

     if (err_weight) cout << "Weights for events are not set. Assumes weights=1 for all" << endl;
     cout << "Writing ROOT file " <<  argv[2] << endl;
     cout << "Total events =" << ntot << endl;
     cout << "Selected events =" << event  << endl;


      return 0;
}
