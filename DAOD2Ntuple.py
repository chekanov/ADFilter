#!/usr/bin/env python
# Convert DAOD to NTuple 
# Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
#
# The goal of this tool is to convert DAOD_PHYS to Ntuple which can be uploaded
# to the ADFilter tool: https://mc.hep.anl.gov/asc/adfilter
# 
# How to run over DAOD_PHYS:
#   lsetup "asetup Athena,main,latest"
#   python ./DAOD2Ntuple.py  --inputlist DAOD_PHYS.37222100._000003.pool.root.1  --outputlist output.root --cmsEnergy 13000
#   
#   S.Chekanov (ANL), Alaettin Serhan Mete (ANL), Wasikul Islam (U.Wicsc)
#  
# The necessary import(s):
import ROOT
import argparse
from array import array
import math

## The main function
def main(filenameinput, filenameoutput, cmsEnergy, cross ):


    MeVtoGeV=0.001

    # typical selection cuts
    PTleptons=30; # PT of leptons 
    PTjets=30;    # PT of jets
    ETAjets=2.4   # Eta of jets
    ETAleptons=2.4 # Eta of leptons

    ifile = ROOT.TFile.Open( filenameinput,  "READ" )
    if not ifile:
        print( "Couldn't open the test input file:", filename )
        return 1

    # counters
    njets=0
    nbjets=0
    nmu=0;
    nel=0;
    nph=0

    # Make a chain of input files
    filelist = ROOT.TChain()
    filelist.AddFile(filenameinput)
    # Read the input file
    evt = ROOT.POOL.TEvent( ROOT.POOL.TEvent.kClassAccess )
    evt.readFrom(filelist) 
    #evt.readFrom(ifile)

    #### ntuple writer ###
    outputFile=ROOT.TFile(filenameoutput, "RECREATE");
    ntuple = ROOT.TTree("Ntuple", "Ntuple for ADFilter");

    N_JET =array('i', [0])
    JET_pt = ROOT.vector('Double32_t')()
    JET_eta = ROOT.vector('Double32_t')()
    JET_phi = ROOT.vector('Double32_t')()
    JET_mass = ROOT.vector('Double32_t')()

    N_bJET = array('i', [0])
    bJET_pt = ROOT.vector('Double32_t')()
    bJET_eta = ROOT.vector('Double32_t')()
    bJET_phi = ROOT.vector('Double32_t')()
    bJET_mass = ROOT.vector('Double32_t')()

    N_EL = array('i', [0])
    EL_pt = ROOT.vector('Double32_t')()
    EL_eta = ROOT.vector('Double32_t')()
    EL_phi = ROOT.vector('Double32_t')()

    N_MU = array('i', [0])
    MU_pt = ROOT.vector('Double32_t')()
    MU_eta = ROOT.vector('Double32_t')()
    MU_phi = ROOT.vector('Double32_t')()

    N_PH = array('i', [0])
    PH_pt = ROOT.vector('Double32_t')()
    PH_eta = ROOT.vector('Double32_t')()
    PH_phi = ROOT.vector('Double32_t')()
    PH_e = ROOT.vector('Double32_t')()

    N_MET = array('i', [0])
    MET_met = ROOT.vector('Double32_t')()
    MET_eta = ROOT.vector('Double32_t')()
    MET_phi = ROOT.vector('Double32_t')()

    Evt_Weight = ROOT.vector('Double32_t')()

    ntuple.Branch( 'JET_n', N_JET, 'N_JET/I' )
    ntuple.Branch( 'JET_pt', JET_pt, 256000, 0 )
    ntuple.Branch( 'JET_eta', JET_eta, 256000, 0 )
    ntuple.Branch( 'JET_phi', JET_phi, 256000, 0 )
    ntuple.Branch( 'JET_mass', JET_mass, 256000, 0 )

    ntuple.Branch( 'N_MET', N_MET, 'N_MET/I' )
    ntuple.Branch( 'MET_met', MET_met, 256000, 0 )
    ntuple.Branch( 'MET_eta', MET_eta, 256000, 0 )
    ntuple.Branch( 'MET_phi', MET_phi, 256000, 0 )

    ntuple.Branch( 'bJET_n', N_bJET, 'N_bJET/I' )
    ntuple.Branch( 'bJET_pt', bJET_pt, 256000, 0 )
    ntuple.Branch( 'bJET_eta', bJET_eta, 256000, 0 )
    ntuple.Branch( 'bJET_phi', bJET_phi, 256000, 0 )
    ntuple.Branch( 'bJET_mass', bJET_mass, 256000, 0 )

    ntuple.Branch( 'MU_n', N_MU, 'N_MU/I' )
    ntuple.Branch( 'MU_pt', MU_pt, 256000, 0 )
    ntuple.Branch( 'MU_eta', MU_eta, 256000, 0 )
    ntuple.Branch( 'MU_phi', MU_phi, 256000, 0 )

    ntuple.Branch( 'PH_n', N_PH, 'N_PH/I' )
    ntuple.Branch( 'PH_pt', PH_pt, 256000, 0 )
    ntuple.Branch( 'PH_eta', PH_eta, 256000, 0 )
    ntuple.Branch( 'PH_phi', PH_phi, 256000, 0 )
    ntuple.Branch( 'PH_e', PH_e, 256000, 0 )

    ntuple.Branch( 'EL_n', N_EL, 'N_EL/I' )
    ntuple.Branch( 'EL_pt', EL_pt, 256000, 0 )
    ntuple.Branch( 'EL_eta', EL_eta, 256000, 0 )
    ntuple.Branch( 'EL_phi', EL_phi, 256000, 0 )

    ntuple.Branch( 'Evt_Weight', Evt_Weight, 256000, 0 )

    meta = ROOT.TH1F("meta","meta",10,0,10);
    meta.Fill("CMS energy [GeV]",cmsEnergy);
    meta.Fill("cross section [PB]",cross);
    meta.SetBinError(1,0)
    meta.SetBinError(2,0)
    meta.Write();

     
    print(f'Total number of events {evt.getEntries()}')
    # https://gist.github.com/kratsg/e5c78c98e666b7d61b86 
    # Loop over the requested number of events
    evtMax = evt.getEntries() #   if flags.Exec.MaxEvents == -1 else flags.Exec.MaxEvents
    for idx in range(0, evtMax):

        JET_pt.clear();
        JET_eta.clear();
        JET_phi.clear();
        JET_mass.clear();

        bJET_pt.clear();
        bJET_eta.clear();
        bJET_phi.clear();
        bJET_mass.clear();

        EL_pt.clear();
        EL_eta.clear();
        EL_phi.clear();

        MU_pt.clear();
        MU_eta.clear();
        MU_phi.clear();

        PH_pt.clear();
        PH_eta.clear();
        PH_phi.clear();
        PH_e.clear();

        MET_eta.clear();
        MET_phi.clear();
        MET_met.clear();

        Evt_Weight.clear();

        if ((idx<100 and idx%10==0) or idx%1000 ==0): print(f'Processing Event #{idx}')
        # Read the current entry
        evt.getEntry(idx)
        # Retrieve the AntiKt4EMPFlowJets jets
        jets = evt.retrieve('xAOD::JetContainer', 'AntiKt4EMPFlowJets')
        # Accessor for the EMFrac of the jet
        jvtAcc = ROOT.SG.ConstAccessor('float')('Jvt')
        # Accessor for the BTagging link of the jet
        btagAcc = ROOT.SG.ConstAccessor('ElementLink<xAOD::BTaggingContainer>')('btaggingLink')
        # Accessor for the GN2v01_pb of the BTag - see how to use this properly instead of auxdataConst below
        #gnAcc = ROOT.SG.ConstAccessor('float')('GN2v01_pb')
        # Loop over the jets to read the pt and EMFrac for demonstration purposes

        # for muons special treatment 
        muonQuality = ROOT.SG.ConstAccessor('unsigned char')('quality')
        muonType = ROOT.SG.ConstAccessor('unsigned short')('muonType')

        #https://ftag.docs.cern.ch/recommendations/algs/r22-preliminary/#gn2v01-b-tagging
        fc=0.2
        ft=0.01
        Ljet,Bjet={},{}
        for jet in jets:
            jvt = None if not jvtAcc.isAvailable(jet) else jvtAcc(jet)

            timing = jet.auxdataConst["float"]("Timing") 
            EMfrac = jet.auxdataConst["float"]("EMFrac")     
            FracSamplingMax = jet.auxdataConst["float"]("FracSamplingMax") 

            pt=MeVtoGeV*jet.pt()
            eta=jet.eta()

            if (EMfrac>0.95): continue  # need to add more here
            if (FracSamplingMax>0.99 and abs(eta)<2): continue; #  high Q factor 
            if (pt<60 and jvt<0.95):                  continue 
            # b-tagging
            btagInfo = None if not btagAcc.isAvailable(jet) else btagAcc(jet)
            score = None if not btagInfo.isValid() else btagInfo.auxdataConst["float"]("GN2v01_pb") # see if this can be done via gnAcc 
            GN2v01_pb= None if not btagInfo.isValid() else btagInfo.auxdataConst["float"]("GN2v01_pb")
            GN2v01_pc= None if not btagInfo.isValid() else btagInfo.auxdataConst["float"]("GN2v01_pc")
            GN2v01_pu= None if not btagInfo.isValid() else btagInfo.auxdataConst["float"]("GN2v01_pu")
            GN2v01_ptau= None if not btagInfo.isValid() else btagInfo.auxdataConst["float"]("GN2v01_ptau")
            DiscriminantB=math.log(GN2v01_pb / (fc*GN2v01_pc + ft*GN2v01_ptau+(1- fc - ft)*GN2v01_pu))

            
            if (pt>PTjets and abs(eta)<ETAjets):
               if (DiscriminantB>0.82): 
                   Bjet[jet]=[pt,eta,jet.phi(),MeVtoGeV*jet.m()] 
                   #print("Bjet 77%",DiscriminantB);
               else:
                   Ljet[jet]=[pt,eta,jet.phi(),MeVtoGeV*jet.m()] 
            #https://ftag.docs.cern.ch/recommendations/algs/r22-preliminary/#gn2v01-b-tagging
            #print(f' >> Jet pt : {jet.pt()} - EMFrac : {emFrac} - GN2v01_pb : {score}')
            #print(f' >> Jet pt : {jet.pt()} - EMFrac {jet.auxdataConst["float"]("EMFrac")}')

        # fill jets
        for i in Ljet.keys():
             kin=Ljet[i]
             JET_pt.push_back(kin[0]);
             JET_eta.push_back(kin[1]);
             JET_phi.push_back(kin[2]);
             JET_mass.push_back(kin[3]);
        for i in Bjet.keys():
             kin=Bjet[i]
             bJET_pt.push_back(kin[0]);
             bJET_eta.push_back(kin[1]);
             bJET_phi.push_back(kin[2]);
             bJET_mass.push_back(kin[3]);
        N_JET[0]=len(JET_pt)
        N_bJET[0]=len(bJET_pt)

        njets=njets+N_JET[0]
        nbjets=nbjets+N_bJET[0]

        electrons = evt.retrieve('xAOD::ElectronContainer', 'Electrons')
        for el in electrons:
            topoetcone40=el.auxdataConst["float"]("topoetcone40")
            tight = ord(el.auxdataConst["char"]("DFCommonElectronsLHTight"))
            if (topoetcone40>0 and tight>0):
                 pt=MeVtoGeV*el.pt()
                 eta=el.eta()
                 if (pt>PTleptons and abs(eta)<ETAleptons):
                    EL_phi.push_back(el.phi())
                    EL_eta.push_back(eta)
                    EL_pt.push_back(pt)
        N_EL[0]=len(EL_pt)
        nel=nel+N_EL[0]


        # Tight, Medium, Loose, VeryLoose, HighPt, LowPt, 
        # corresponding to the "MuQuality" property set to 0, 1, 2, 3, 4 and 5
        # https://twiki.cern.ch/twiki/bin/view/Atlas/MuonSelectionTool
        xAOD_Muon_Combined=1 # Combined 
        xAOD_MuonQuality =3  # Medium ??  
        muons = evt.retrieve('xAOD::MuonContainer', 'Muons')
        for mu in muons:
            mutype =int(None if not muonType.isAvailable(mu) else muonType(mu))
            quality = ord(None if not muonQuality.isAvailable(mu) else muonQuality(mu))
            #print("Type=",mutype,"Quality=",quality)
            if (quality == xAOD_MuonQuality and mutype == xAOD_Muon_Combined):
               pt=MeVtoGeV*mu.pt()
               eta=mu.eta()
               if (pt>PTleptons and abs(eta)<ETAleptons):
                 MU_phi.push_back(mu.phi())
                 MU_eta.push_back(eta)
                 MU_pt.push_back(pt)
        N_MU[0]=len(MU_pt)
        nmu=nmu+N_MU[0]

        photons = evt.retrieve('xAOD::PhotonContainer', 'Photons')
        for ph in photons:
            topoetcone40=1 # h.auxdataConst["float"]("topoetcone40")
            tight=ord(ph.auxdataConst["char"]("Tight"))
            if (topoetcone40>0 and tight>0):
                    pt=MeVtoGeV*ph.pt()
                    eta=ph.eta()
                    if (pt>PTleptons and abs(eta)<ETAleptons): 
                      PH_phi.push_back(ph.phi())
                      PH_eta.push_back(eta)
                      PH_pt.push_back(pt)
        N_PH[0]=len(PH_pt)
        nph=nph+N_PH[0]

        met = evt.retrieve('xAOD::MissingETContainer', 'MET_Core_AntiKt4EMPFlow')
        for m in met:
              x=MeVtoGeV*m.mpx()
              y=MeVtoGeV*m.mpy()
              phi = math.atan2(y, x)
              pt=math.sqrt(x*x+y*y)
              MET_phi.push_back(phi )
              MET_eta.push_back(0)
              MET_met.push_back(pt)
        N_MET[0]=len(MET_met)
    
        # weights 
        Evt_Weight.push_back(1)
        Evt_Weight.push_back(1)

        # fill ntuple
        ntuple.Fill()

    #ifile.Close()
    print("Write the TTree to the output file:",filenameoutput);
    #outputFile.Write("",TFile.kOverwrite)
    #ntuple.cd()
    #ntuple.SetDirectory(outputFile);
    outputFile.WriteObject(ntuple,"Ntuple")
    #ntuple.Write()
    #ntuple.Show(1)
    #Close the output file
    outputFile.Close();
    print("",njets)
    print("Output summary:",njets)
    print("Nr of light jets=",njets)
    print("Nr of b-tag jets=",nbjets)
    print("Nr of electrons =",nel)
    print("Nr of muons =",nmu)
    print("Nr of photons =",nph)
 

# Run the main() function:
if __name__ == "__main__":
   import sys

   parser = argparse.ArgumentParser(
       description="Extracts a few basic quantities from the xAOD file and dumps them into a text file")
   parser.add_argument("--outputlist", help="List of output ROOT files",
                       nargs='+', action="store", default=False)
   parser.add_argument("--cmsEnergy", help="CMS energy in GeV",
                       nargs='+', action="store", default=False)
   parser.add_argument("--inputlist", help="List of  DAOD_PHYS files",
                       nargs='+', action="store", default=False)
   parser.add_argument("--crossSectionPB", help="Cross section in [pb]",
                       nargs='+', action="store", default=False)

   args = parser.parse_args()

   cmsEnergy=13000
   filename="output.root"
   cross=1.0;
   if args.outputlist:
            filename=args.outputlist[0]
   if args.inputlist:
            filenameinput=args.inputlist[0]
   if args.cmsEnergy:
            print("CMS energy [GeV] =",args.cmsEnergy[0]);
            cmsEnergy=float(args.cmsEnergy[0])
   if args.crossSectionPB:
            print("Cross section [pb] =",args.crossSectionPB[0]);
            cross=float(args.crossSectionPB[0])

   sys.exit( main( filenameinput, filename, cmsEnergy, cross  ) )


