# Anomaly detection https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.10.3542&rep=rep1&type=pdf
# Task: Process events using the model
# S.Chekanov (ANL)

import os,sys
import random
####*IMPORANT*: Have to do this line *before* importing tensorflow
os.environ['PYTHONHASHSEED']=str(1)

# Is to create a tree too?
save_tree=False 

sys.path.append("modules/")
from global_module import *
from ROOT import TH1D, TF1, TProfile2D, TEllipse, THStack, TRandom3, TFile, TTree, TLatex, TLegend, TPaveText, TGraphErrors, kRed, kBlue, kGreen, kCyan, kAzure, kYellow, kTRUE
import ROOT

import numpy
import pandas
import matplotlib
import seaborn
import tensorflow
import pickle
print('Numpy version      :' , numpy.__version__)
print('Pandas version     :' ,pandas.__version__)
print('Matplotlib version :' ,matplotlib.__version__)
print('Seaborn version    :' , seaborn.__version__)
print('Tensorflow version :' , tensorflow.__version__)


import numpy as np
import pandas as pd
#pd.set_option('display.max_columns', None)
#pd.set_option('display.max_row', None)
#import matplotlib.pyplot as plt
#plt.rcdefaults()
#from pylab import rcParams
#import seaborn as sns
import datetime
# import matplotlib
# matplotlib.use('Agg') # set the backend before importing pyplo. Fix Invalid DISPLAY variable 
# from matplotlib import pyplot as plt
####### Deep learning libraries
import tensorflow as tf
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, Dense


# Data Preprocessing
import pandas
import matplotlib
import seaborn
import tensorflow
import tensorflow as tf

# random seeds fixed
RANDOM_SEED = 101
os.environ['PYTHONHASHSEED']=str(1)
tf.random.set_seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED*2)
random.seed(RANDOM_SEED*3)
print("Use fixed seed=",RANDOM_SEED)
os.environ["OMP_NUM_THREADS"] = "1"
physical_devices = tf.config.list_physical_devices('CPU')
tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(1)


print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
n = len(sys.argv)
if (n != 6):
      print ("No arguments!. Need at least 5: model, events, input, output, working_point") 
      sys.exit()

model=sys.argv[1]
MaxEvents=int(sys.argv[2]) 
inputData=sys.argv[3]
outputData=sys.argv[4]

workingPoint=float(sys.argv[5]) 

# get working points and teh cut on the loss..
CutOutlierMC=CutOutlier_10PB 

KSoutfiles=inputData.replace(".root","");

rootfile=KSoutfiles+"_ADFilter.root"

print("-> Run over max events=",MaxEvents)
print("-> Use model =",model)
print("-> Input data =",inputData)
print("-> Output  data =",outputData)
print("-> Working point =",workingPoint,"Corresponds to the cut on the loss=",CutOutlierMC)
print("-> Output ROOT file=",rootfile)
print("")

# if data and only 10%
onlyFraction=1.0
if (outputData.find("1percent")>-1):
                     onlyFraction=0.01
                     print("Process fraction = ", onlyFraction*100,"%")
if (outputData.find("10percent")>-1):
                     onlyFraction=0.1
                     print("Process fraction =", onlyFraction*100,"%")


inRMM=[inputData]
inLabel=[outputData]


proc=[ inputData ]
outRMM=[ outputData ]

rfile=[]
for i in proc:
     rfile.append(ROOT.TFile.Open(i))
     print(i)


# double lit to keep data for dataframe
columnsX=[]
for i in range(1,mSize*mSize+1):
       columnsX.append("V_"+str(i))
# last column labels the data (put 0) 
# columnsX.append("Label")
df = pd.DataFrame(columns=columnsX)
print("DF size=",df.size," DF shape=",df.shape," DF dimension=",df.ndim)



"""
# load json and create model
fj1="figs/model.json"
print("--> Read = ",fj1) 
json_file = open(fj1, 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json, custom_objects={'leaky_relu': tf.nn.leaky_relu})
# load weights into new model
fj2="figs/model.h5"
print("--> Read = ",fj2)
loaded_model.load_weights(fj2)
print("--> Loaded model from disk")
"""


fj1="models"
loaded_model = tf.keras.models.load_model( fj1 )
print("--> Loaded model from "+fj1)


loaded_model.summary()
loaded_model.compile(optimizer='adam', loss='mse' )


# apply Standardization and MinMax?
IsStandard=False

if (IsStandard):
  print("")
  print("Data Standardization.. so that the mean of observed values is 0 and the standard deviation is 1.");
  import pickle
  scaler_filename = fj1+"/StandardScaler.pkl"
  print("Read StandardScaler =",scaler_filename)
  scalerStandard = pickle.load(open(scaler_filename, 'rb'))
  print ("Data scaling.. Can be skipped since RMM [0-1]. But you ran standartisation before!")
  # Data Scaling
  scaler_filename = fj1+"/MinMaxScaler.pkl"
  print("Read fitted MinMaxScaler =",scaler_filename)
  scalerMinMax = pickle.load(open(scaler_filename, 'rb'))
  print("")


IsReadCommonEmptyColumns=1
# 1 drop columns based on common vector
# 2 drop columns as found by the current dataframe
file0=""
if (IsReadCommonEmptyColumns==1):
   file0="columns_with_0.txt"
if (IsReadCommonEmptyColumns==2):
   file0=fj1+"columns_with_0.txt"
print ("Read columns with 0 from ",file0)
dcol0=pd.read_csv(file0,header = None)
col0=dcol0[dcol0.columns[0]]

ka="";
# z-score 
h_loss=TH1D("Loss"+ka,"Loss distribution ( ln(Loss) ) "+ka,200,-13, -3)
h_debug=TH1D("Eventflow"+ka,"Event flow"+ka,10,0,10)
h_cut=TH1D("LossCut"+ka,"Cut on the ln( Loss )"+ka,200,-13, -3)

binsM = TH1D("bins_m"+ka, "bins_m"+ka, len(mjjBins)-1, mjjBins);
for j in range( len(mjjBins)-1):
       x=mjjBins[j+1]-mjjBins[j];
       binsM.Fill(mjjBins[j]+0.5*x,x);

# masses 
h1=TH1D("Mjj"+ka,"jet-jet mass  in AR"+ka, len(mjjBins)-1, mjjBins )
h2=TH1D("Mbb"+ka,"bjet - bjet mass in AR"+ka, len(mjjBins)-1, mjjBins)
h3=TH1D("Mjb"+ka,"jet - bjet mass in AR"+ka, len(mjjBins)-1, mjjBins)
h4=TH1D("Mee"+ka,"ee mass in AR"+ka,200,1,401)
h5=TH1D("Mmm"+ka,"muon-muon mass in AR "+ka,200,1,401)
h6=TH1D("Mje"+ka,"jet-e mass in AR"+ka, len(mjjBins)-1, mjjBins)
h7=TH1D("Mjm"+ka,"jet-muon mass in AR"+ka, len(mjjBins)-1, mjjBins)
h8=TH1D("Mjg"+ka,"jet-gamma mass in AR"+ka, len(mjjBins)-1, mjjBins)
h9=TH1D("Mgg"+ka,"gamma-gamma mass in AR"+ka, len(mjjBins)-1, mjjBins)
h10=TH1D("Mge"+ka,"gamma-e mass in AR"+ka,len(mjjBins)-1, mjjBins)
h11=TH1D("Mbg"+ka,"bjet-gamma mass in AR"+ka,len(mjjBins)-1, mjjBins)
h12=TH1D("Mgm"+ka,"gamma-muon mass in AR"+ka,len(mjjBins)-1, mjjBins)
h13=TH1D("Mem"+ka,"e-muon mass in AR"+ka,100,1,1001)
h14=TH1D("Mbe"+ka,"bjet-e mass in AR"+ka,len(mjjBins)-1, mjjBins)
h15=TH1D("Mbm"+ka,"bjet-muon in AR"+ka,len(mjjBins)-1, mjjBins)


# before NN cut 
h1b=TH1D("Mjj_b"+ka,"jet-jet before Autoencoder"+ka, len(mjjBins)-1, mjjBins )
h2b=TH1D("Mbb_b"+ka,"bjet - bjet before Autoencoder"+ka, len(mjjBins)-1, mjjBins)
h3b=TH1D("Mjb_b"+ka,"jet - bjet before Autoencoder"+ka, len(mjjBins)-1, mjjBins)
h4b=TH1D("Mee_b"+ka,"ee mass before Autoencoder"+ka,200,1,401)
h5b=TH1D("Mmm_b"+ka,"muon-muon before Autoencoder"+ka,200,1,401)
h6b=TH1D("Mje_b"+ka,"jet-e before Autoencoder"+ka, len(mjjBins)-1, mjjBins)
h7b=TH1D("Mjm_b"+ka,"jet-muon before Autoencoder"+ka, len(mjjBins)-1, mjjBins)
h8b=TH1D("Mjg_b"+ka,"jet-gamma before Autoencoder"+ka, len(mjjBins)-1, mjjBins)
h9b=TH1D("Mgg_b"+ka,"gamma-gamma mass before Autoencoder"+ka, len(mjjBins)-1, mjjBins)
h10b=TH1D("Mge_b"+ka,"gamma-e mass before Autoencoder"+ka,len(mjjBins)-1, mjjBins)
h11b=TH1D("Mbg_b"+ka,"bjet-gamma before Autoencoder"+ka,len(mjjBins)-1, mjjBins)
h12b=TH1D("Mgm_b"+ka,"gamma-muon mass before Autoencoder"+ka,len(mjjBins)-1, mjjBins)
h13b=TH1D("Mem_b"+ka,"e-muon mass before Autoencoder"+ka,100,1,1001)
h14b=TH1D("Mbe_b"+ka,"bjet-e mass before Autoencoder"+ka,len(mjjBins)-1, mjjBins)
h15b=TH1D("Mbm_b"+ka,"bjet-muon mass before Autoencoder"+ka,len(mjjBins)-1, mjjBins)



if save_tree:
    tree = ROOT.TTree("output","output")
    process = array('i', [0])
    mc_sf = array('f', [0.])
    Mjj = array('f', [0.])
    Mbb = array('f', [0.])
    Mjb = array('f', [0.])
    Mee = array('f', [0.])
    Mmm = array('f', [0.])
    Mje = array('f', [0.])
    Mjm = array('f', [0.])
    Mjg = array('f', [0.])
    Mbe = array('f', [0.])
    Mbm = array('f', [0.])
    Mbg = array('f', [0.])
    tree.Branch("process", process, 'process/I')
    tree.Branch("mc_sf", mc_sf, 'mc_sf/F')
    tree.Branch("Mjj", Mjj, 'Mjj/F')
    tree.Branch("Mbb", Mbb, 'Mbb/F')
    tree.Branch("Mjb", Mjb, 'Mjb/F')
    tree.Branch("Mee", Mee, 'Mee/F')
    tree.Branch("Mmm", Mmm, 'Mmm/F')
    tree.Branch("Mje", Mje, 'Mje/F')
    tree.Branch("Mjm", Mjm, 'Mjm/F')
    tree.Branch("Mjg", Mjg, 'Mjg/F')
    tree.Branch("Mbe", Mbe, 'Mbe/F')
    tree.Branch("Mbm", Mbm, 'Mbm/F')
    tree.Branch("Mbg", Mbg, 'Mbg/F')




h1.Sumw2();
h2.Sumw2();
h3.Sumw2();
h4.Sumw2();
h5.Sumw2();
h6.Sumw2();
h7.Sumw2();
h8.Sumw2();
h9.Sumw2();
h10.Sumw2();
h11.Sumw2();
h12.Sumw2();
h13.Sumw2();
h14.Sumw2();
h15.Sumw2();

h1b.Sumw2();
h2b.Sumw2();
h3b.Sumw2();
h4b.Sumw2();
h5b.Sumw2();
h6b.Sumw2();
h7b.Sumw2();
h8b.Sumw2();
h9b.Sumw2();
h10b.Sumw2();
h11b.Sumw2();
h12b.Sumw2();
h13b.Sumw2();
h14b.Sumw2();
h15b.Sumw2();

# max number of masses to be analysed
MaxNumberOfMasses=11


ntot=0
kk=0;
events=0;
inputs=0
outputs=0

mean={}
sigma={}
events={}


xfsum1=inRMM[0].replace(".zip","_summary.txt")
print ("Read summary file: pos, tot,av,sigma=",xfsum1) 

# signal is avaluated with respect data (or SM MC) 
if (ka.find("signal")>-1):
         xfsum1="data/data14invfb_summary.txt"
         print ("Read summary file: pos, tot,av,sigma=",xfsum1) 

ntot=0
print ("Start processing..") 
# how many chunks with RMM 
evtINchunk=1000
nchunk=0
chunk=0
evt=0

for i in range(len(proc)):
   ev=0

   # initialize
   RMM = np.zeros(shape=(evtINchunk, mSize*mSize))
   masses=np.zeros(shape=(evtINchunk, MaxNumberOfMasses))
   weights=[]
   runs=[]
   events=[]

   # total events in file
   NtotInFile=(rfile[i].inputNN).GetEntries()
   print("Analyse ",NtotInFile," from file=",rfile[i])

   for event in rfile[i].inputNN:

       if (onlyFraction<1.0):
          ran=random.uniform(0, 1)
          if (ran>onlyFraction): continue

       NN=(event.proj).size()
       a=event.proj
       inx1=event.proj_index1
       inx2=event.proj_index2
       Trun = event.run
       Tevent=event.event
       Tweight=event.weight # for MC with weigths
       weight=Tweight;

       #ST='%.5E' % Decimal(Tweight)
       #pos=str(Trun)+"#"+str(Tevent)+"#"+str(ST)

       emptyMatrix = numpy.zeros(shape=(mSize,mSize))
       txt=""

       v_mjj  = 0;
       v_mbb  = 0;
       v_mjb  = 0;
       v_mee  = 0;
       v_mmm  = 0;
       # additional
       v_mje  = 0;
       v_mjm  = 0;
       v_mjg  = 0;
       v_mbe  = 0;
       v_mbm  = 0;
       v_mbg  = 0;


       for i3 in range(NN):
              w=inx1[i3];
              h=inx2[i3];
              val=float(a[i3])
              emptyMatrix[w][h] = val
              # mjj
              if (h==mjj[0] and w==mjj[1]): v_mjj=val*CMS
              # mbb
              if (h==mbb[0] and w==mbb[1]): v_mbb=val*CMS
              # mee
              if (h==mee[0] and w==mee[1]): v_mee=val*CMS
              #mumu
              if (h==mmumu[0] and w==mmumu[1]): v_mmm=val*CMS
              # mjb 
              if (h==mbj[0] and w==mbj[1]):  v_mjb=val*CMS
              # mje 
              if (h==mje[0] and w==mje[1]):  v_mje=val*CMS
              # mje 
              if (h==mjmu[0] and w==mjmu[1]):  v_mjm=val*CMS
              # mjg 
              if (h==mjg[0] and w==mjg[1]):  v_mjg=val*CMS
              # mbe 
              if (h==mbe[0] and w==mbe[1]):  v_mbe=val*CMS
              # mje 
              if (h==mbmu[0] and w==mbmu[1]):  v_mbm=val*CMS
              # mjg 
              if (h==mbg[0] and w==mbg[1]):  v_mbg=val*CMS


       # fill masses before any cut
       h1b.Fill( v_mjj, weight )
       h2b.Fill( v_mbb, weight )
       h3b.Fill( v_mjb, weight )
       h4b.Fill( v_mee, weight )
       h5b.Fill( v_mmm, weight )
       # additional
       h6b.Fill( v_mje, weight ) # j+e 
       h7b.Fill( v_mjm, weight ) # j+mu 
       h8b.Fill( v_mjg, weight ) # j+gamma 
       h14b.Fill( v_mbe, weight ) # b+e 
       h15b.Fill( v_mbm, weight ) # b+mu 
       h11b.Fill( v_mbg, weight ) # b+gamma 


       # print(v_mje,  v_mjm) 

       # flatten
       dataRMM=(emptyMatrix.flatten()).tolist()
       RMM[evt,:]=dataRMM       
       masses[evt,:]=numpy.array([v_mjj,v_mbb,v_mjb,v_mee,v_mmm,v_mje,v_mjm,v_mjg,v_mbe,v_mbm,v_mbg])
       weights.append(weight) 
       runs.append(Trun)
       events.append(Tevent)

       evt=evt+1   # events in chunk 
       ev=ev+1     # events in this file 
       ntot=ntot+1 # all events 
       if (ev ==  NtotInFile or evt%evtINchunk==0):

                     df = pd.DataFrame(data=RMM, columns=columnsX)
                     df=df.drop(col0, axis = 1)
                     #print("Apply scalers and remove 0 columns: DF size=",df.size," DF shape=",df.shape," DF dimension=",df.ndim)
                     RMM_T = df.to_numpy()

                     if (IsStandard):
                       RMM_T = scalerStandard.transform(RMM_T)
                       RMM_T = scalerMinMax.transform(RMM_T)

                     predictions = loaded_model.predict( RMM_T )
                     train_loss = tf.keras.losses.mae(predictions, RMM_T).numpy()

                     nle=len(train_loss)
                     if (ev ==  NtotInFile):
                                          nle=ev- nchunk*evtINchunk 
                                          print(" -> Last event ",ev, " from ",NtotInFile, "remaining=",nle)

                     for ch in range(nle):

                            xloss=train_loss[ch]
                            we=weights[ch]
                            xrun=runs[ch]
                            xevent=events[ch]

                            xlog= math.log(xloss) 
                            #print(xlog)
                            h_loss.Fill(xlog, we)
                            h_debug.Fill("Input events",1)
                            h_cut.Fill(CutOutlierMC);

                            if (xlog < CutOutlierMC): continue # reject SM using MC outlier 
                            # This is side-band control region -10.10 - 10.0
                            if (workingPoint==-1):
                                               if (xlog >  CutOutlierMC+0.10): continue  

                            h_debug.Fill("Events in AR",1)
                            ma=masses[ch]
                            # after SM reject
                            h1.Fill( ma[0], we )
                            h2.Fill( ma[1], we )
                            h3.Fill( ma[2], we )
                            h4.Fill( ma[3], we )
                            h5.Fill( ma[4], we )
                            # additional masses
                            h6.Fill( ma[5], we )
                            h7.Fill( ma[6], we  )
                            h8.Fill( ma[7], we )
                            h14.Fill( ma[8], we )
                            h15.Fill( ma[9], we )
                            h11.Fill(ma[10], we )

                            if save_tree:
                                process[0] = i
                                mc_sf[0] = we
                                Mjj[0] = ma[0]
                                Mbb[0] = ma[1]
                                Mjb[0] = ma[2]
                                Mee[0] = ma[3]
                                Mmm[0] = ma[4]
                                Mje[0] = ma[5]
                                Mjm[0] = ma[6]
                                Mjg[0] = ma[7]
                                Mbe[0] = ma[8]
                                Mbm[0] = ma[9]
                                Mbg[0] = ma[10]
                                tree.Fill()



                            #print(ma[5], ma[6]) 

                     if (MaxEvents>0):
                          if (ntot>MaxEvents):
                             print ("Stop loop ",MaxEvents); break;

                     # reset
                     RMM = np.zeros(shape=(evtINchunk, mSize*mSize))
                     masses=np.zeros(shape=(evtINchunk, MaxNumberOfMasses))
                     weights=[]
                     runs=[]
                     events=[]
                     print("Fill chunk ",nchunk," with ",evtINchunk, " events. Tot=",ev,"  mjj=",int(h1.GetEntries()))
                     nchunk=nchunk+1
                     evt=0
       if (MaxEvents>0):
                 if (ntot>MaxEvents):
                   print ("Finish after ",MaxEvents, " events"); break;
 

print("Total events =",ev);


import os.path
if (os.path.exists(rootfile) == False):
  print("File=",rootfile," does not exist. Make it")
  hfile=TFile(rootfile,"RECREATE","signatures")
  h_debug.Write()
  hfile.Close()


hfile=TFile(rootfile,"RECREATE","signatures")
h_debug.Write()
h_loss.Write()
h_cut.Write()
h1.Write()
h2.Write()
h3.Write()
h4.Write()
h5.Write()
h6.Write()
h7.Write()
h8.Write()
h9.Write()
h10.Write()
h11.Write()
h12.Write()
h13.Write()
h14.Write()
h15.Write()
# before
h1b.Write()
h2b.Write()
h3b.Write()
h4b.Write()
h5b.Write()
h6b.Write()
h7b.Write()
h8b.Write()
h9b.Write()
h10b.Write()
h11b.Write()
h12b.Write()
h13b.Write()
h14b.Write()
h15b.Write()
# bins
binsM.Write()

# tree
if save_tree:
    tree.Write()

hfile.Print()

hfile.Close()
print ("Write=",rootfile) 





