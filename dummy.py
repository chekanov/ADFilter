#!/usr/bin/env python
#
# Create Dummy record for ADFilter
# python ./dummy.py --outputlist output.root  --cmsEnergy 13000 
# You should have ROOT with python compiled.
# S.V.Chekanov (ANL) 

# The necessary import(s):
import ROOT
import argparse
from array import array


## C/C++ style main function
def main( filename, cmsENERGY ):

    # The name of the application:
    APP_NAME = "dummy"
    # number of events
    NEvents=1000

    # Set up a logger object:
    import logging
    logger = logging.getLogger( APP_NAME )
    logger.setLevel( logging.INFO )
    hdlr = logging.StreamHandler( sys.stdout )
    frmt = logging.Formatter( "%(name)-14s%(levelname)8s %(message)s" )
    hdlr.setFormatter( frmt )
    logger.addHandler( hdlr )

    #### ntuple writer ###
    cmsEnergy=cmsENERGY;
    outputFileName=filename
    outputFile=ROOT.TFile(outputFileName, "RECREATE");
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
    meta.Write();

    print("Total number of events=", NEvents)
    for entry in range( NEvents ):
        if (entry%10 == 0): logger.info( "Processing entry %i" % entry )
        # Print the properties of the electrons:


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

        # define electrons
        N_EL[0] = 1;
        for j in  range(N_EL[0]):
            EL_phi.push_back(0.0)
            EL_eta.push_back(0.0) 
            EL_pt.push_back(100.0) 

        # define muons
        N_MU[0] =  1 
        for j in range(N_MU[0]):
            MU_phi.push_back(0.0)
            MU_eta.push_back(0.0)
            MU_pt.push_back(100.0)

        # define your photons
        N_PH[0] = 1 
        for j in range(N_PH[0]):  
            PH_phi.push_back(0.0)
            PH_eta.push_back(0.0)
            PH_pt.push_back(100.0)
            PH_e.push_back(0.0)

        # define here light-flavour jets
        N_JET[0] = 1
        for j in range(N_JET[0]):
            JET_phi.push_back(0)
            JET_eta.push_back(0.0)
            JET_pt.push_back(100.0)
            JET_mass.push_back(100.0)

        # define here b-jets
        N_bJET[0] = 1 
        for j in range(N_bJET[0]):
            bJET_phi.push_back(0.0)
            bJET_eta.push_back(0.0)
            bJET_pt.push_back(100.0)
            bJET_mass.push_back(100.0)

        # define MET. 
        N_MET[0] = 1 # should be 1 always! 
        for j in range(N_MET[0]):
            MET_phi.push_back(0.0)
            MET_eta.push_back(0.0)
            MET_met.push_back(0.0)

        # define event weight
        Evt_Weight.push_back(1)
        Evt_Weight.push_back(1)

        ntuple.Fill()
           


    print("Write the TTree to the output file:",outputFileName);
    #outputFile.Write("",TFile.kOverwrite)
    ntuple.Write()
    #ntuple.Show(1)

    # Close the output file
    outputFile.Close();
    print("Write=",outputFileName)


    # Return gracefully:
    return 0

# Run the main() function:
if __name__ == "__main__":
   import sys

   parser = argparse.ArgumentParser(
       description="Extracts a few basic quantities from the xAOD file and dumps them into a text file")
   parser.add_argument("--outputlist", help="Optional list of output ROOT files",
                       nargs='+', action="store", default=False)
   parser.add_argument("--cmsEnergy", help="Optional CMS energy",
                       nargs='+', action="store", default=False)

   args = parser.parse_args()

   cmsEnergy=13000
   filename="dummy.root"
   if args.outputlist:
            filename=args.outputlist[0] 
   if args.cmsEnergy:
            print("CMS energy =",args.cmsEnergy[0]);
            cmsEnergy=float(args.cmsEnergy[0]) 

   sys.exit( main( filename, cmsEnergy  ) )

