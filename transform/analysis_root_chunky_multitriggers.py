# Anomaly detection https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.10.3542&rep=rep1&type=pdf
# Task: Process events using the model
# allows to use AE, PAE and HAE
# S.Chekanov (ANL)

import os,sys
import random
####*IMPORANT*: Have to do this line *before* importing tensorflow
os.environ['PYTHONHASHSEED']=str(1)

SCRIPT_DIR=os.getcwd()
print("Script dirrectory=",SCRIPT_DIR)


# Is to create a tree too?
save_tree=False 

#sys.path.append("modules/")
#from AtlasStyle import *
#from AtlasUtils import *
#from global_module import *
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


# Anomaly detection https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.10.3542&rep=rep1&type=pdf

print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
n = len(sys.argv)
if (n != 8):
      print ("No arguments!. Need at least 8 parameters: model, events, input, output, type, nn_type working_point") 
      sys.exit()

model=sys.argv[1]
MaxEvents=int(sys.argv[2]) 
inputData=sys.argv[3]
outputData=sys.argv[4]
# trigger types
TYPE=sys.argv[5]

#  AE type
NN_TYPE=sys.argv[6]

workingPoint=-1
if sys.argv[7] != "SB":
    workingPoint=float(sys.argv[7]) 

#
# load anomaly region AR
ARfile="AR_AE.txt"
if (NN_TYPE == "HAE"): ARfile=SCRIPT_DIR+"/models_multitriggers/"+"AR_HAE.txt"
if (NN_TYPE == "HAE_RUN2"): ARfile=SCRIPT_DIR+"/models_multitriggers/"+"AR_HAE_RUN2.txt"
if (NN_TYPE == "HAE_RUN3"): ARfile=SCRIPT_DIR+"/models_multitriggers/"+"AR_HAE_RUN3.txt"

# AE trained on mix of Run2+Run3 data
if (NN_TYPE == "HAE_RUN23"): ARfile=SCRIPT_DIR+"/models_multitriggers/"+"AR_HAE_RUN23.txt"
# Run3 using Run23 AE
if (NN_TYPE == "HAE_RUN23_4RUN3"): ARfile=SCRIPT_DIR+"/models_multitriggers/"+"AR_HAE_RUN23.txt"
# Run2 using Run23 AE
if (NN_TYPE == "HAE_RUN23_4RUN2"): ARfile=SCRIPT_DIR+"/models_multitriggers/"+"AR_HAE_RUN23.txt"

print("Using AR from = ",ARfile)

import json
import os.path
ARvalues={}
if  os.path.isfile(ARfile):
    ARvalues = json.load(open( ARfile ))
    print("Anomaly regions from the file=",ARfile)
    print(ARvalues)
else:
    print("Exit. We cannot read AR file from ", ARfile);
    sys.exit()


# this is input parameter (integer) that tells which data stream
input_trigger_type=int(TYPE.replace("t",""));

# when trigger type is complex t1_t2. This means t1 is used as input, but autoencoder is taken from t2
# this is for cross validation, when input for NN was taken from different stream
OUT_TRIGGER=TYPE
if (TYPE.find("_")>-1):
              datastream=TYPE.split("_")
              TYPE=datastream[0]        #  type for input data 
              APPLIED_TYPE=datastream[1] # applied to this stream 
              input_trigger_type=int(APPLIED_TYPE.replace("t",""));
              print("Cross validation! Trigger type for AE = ", input_trigger_type)
              print("Cross validation! Input data from  = ", TYPE)

# get working points and then cut on the loss..
# When we use 10 pb working point, we assume we want to use Data AR
# get working points and teh cut on the loss..
CutOutlierMC=CutOutlier_10PB


KSoutfiles=inputData.replace(".root","");


rootfile=KSoutfiles+"_ADFilter.root"


# move to MC 
if (outputData.find("mc20")>-1): rootfile=rootfile.replace("root/"+outputData,"root/mc20/"+outputData); 
# new name for Run3
if (inputData.find("data2022")>-1): rootfile=rootfile.replace("percent","percent_run3");
# new name for Run3
if (inputData.find("data2023")>-1): rootfile=rootfile.replace("percent","percent_run3");



print("-> Run over max events=",MaxEvents)
print("-> Use model =",model)
print("-> Input data =",inputData)
print("-> Output  data =",outputData)
print("-> Trigger type for INPUT =",TYPE)
print("-> Trigger type for Autoencoder =",input_trigger_type)
print("-> NN   type =",NN_TYPE)
print("-> Working point =",workingPoint,"Corresponds to the cut on the loss=",CutOutlierMC)
print("-> Output ROOT file=",rootfile)
print("")

# if you want to test training against next (different!) trigger stream
# it should show anomaly 
# input_trigger_type=(input_trigger_type%7)+1

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

def custom_mse_loss(y_true, y_pred):
    # Assuming y_true and y_pred are tensors of shape (batch_size, N)
    n = tf.shape(y_true)[-1]  # N, original dimension of the input
    # Select N-1 dimensions for computation of MSE
    reduced_dim = n - 1
    y_true_reduced = y_true[:, :reduced_dim]
    y_pred_reduced = y_pred[:, :reduced_dim]
    mse_loss = tf.keras.losses.mse(y_true_reduced, y_pred_reduced)
    return mse_loss


loaded_model = None 
if (NN_TYPE == "AE"):
  fj1="models/training_v11/"+TYPE+"/"+model+"/mdAEleakyRelu800_400_200_400_800_bs100_sc1.0_dbFalse_seed101_laten20_lr0.001_sergei/models"
  print("--> Loading AE model from "+fj1)
  loaded_model = tf.keras.models.load_model( fj1 )
  loaded_model.summary()
  loaded_model.compile(optimizer='adam', loss='mse' )

if (NN_TYPE == "PAE"):
  fj1="models/training_v11/parametrised/"+model+"/mdAEleakyRelu800_400_200_400_800_bs100_sc1.0_dbFalse_seed101_laten20_lr0.001_sergei/models"
  print("--> Loading  PARAMETERIZED model from "+fj1)
  loaded_model = tf.keras.models.load_model(fj1, custom_objects={'custom_mse_loss': custom_mse_loss})
  loaded_model.summary()
  loaded_model.compile(optimizer='adam', loss='mse' )
# one-hot-encoder
if (NN_TYPE == "HAE"):
  fj1="models/training_v11/one_hot_encoding/"+model+"/mdAEleakyRelu800_400_200_400_800_bs100_sc1.0_dbFalse_seed101_laten20_lr0.001_sergei/models"
  print("--> Loading ONE-HOT-ENCORDER  model from "+fj1)
  loaded_model = tf.keras.models.load_model(fj1, custom_objects={'custom_mse_loss': custom_mse_loss})
  loaded_model.summary()
  loaded_model.compile(optimizer='adam', loss='mse' )
if (NN_TYPE == "HAE_RUN2"):
  fj1="models/training_v12/run2/one_hot_encoding/"+model+"/mdAEleakyRelu800_400_100_400_800_bs100_sc1.0_dbFalse_seed101_laten20_lr0.001_sergei/models"
  print("--> Loading ONE-HOT-ENCORDER  for RUN2 model from "+fj1)
  loaded_model = tf.keras.models.load_model(fj1, custom_objects={'custom_mse_loss': custom_mse_loss})
  loaded_model.summary()
  loaded_model.compile(optimizer='adam', loss='mse' )
if (NN_TYPE == "HAE_RUN3"):
  fj1=SCRIPT_DIR+"/models_multitriggers/training_v12/run3/one_hot_encoding/"+model+"/mdAEleakyRelu800_400_100_400_800_bs100_sc1.0_dbFalse_seed101_laten20_lr0.001_sergei/models"
  print("--> Loading ONE-HOT-ENCORDER  for RUN3 model from "+fj1)
  loaded_model = tf.keras.models.load_model(fj1, custom_objects={'custom_mse_loss': custom_mse_loss})
  loaded_model.summary()
  loaded_model.compile(optimizer='adam', loss='mse' )
if (NN_TYPE.find("HAE_RUN23")>-1): # training using Run2+Run3 
  fj1=SCRIPT_DIR+"/models_multitriggers/training_v12/run23/one_hot_encoding/"+model+"/mdAEleakyRelu800_400_100_400_800_bs100_sc1.0_dbFalse_seed101_laten20_lr0.001_sergei/models"
  print("--> Loading ONE-HOT-ENCORDER  for RUN2+RUN3 model from "+fj1)
  loaded_model = tf.keras.models.load_model(fj1, custom_objects={'custom_mse_loss': custom_mse_loss})
  loaded_model.summary()
  loaded_model.compile(optimizer='adam', loss='mse' )

if (loaded_model == None):
              print("Model was not loaded!")
              sys.exit()

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
   file0=SCRIPT_DIR+"/columns_with_0.txt"
if (IsReadCommonEmptyColumns==2):
   file0=fj1+"columns_with_0.txt"
print ("Read columns with 0 from ",file0)
dcol0=pd.read_csv(file0,header = None)
col0=dcol0[dcol0.columns[0]]


ka=inLabel[0]
# z-score 
h_loss=TH1D("Loss_"+ka,"Loss_"+ka,200,-13, -3)

binsM = TH1D("bins_m_"+ka, "bins_m_"+ka, len(mjjBins)-1, mjjBins);
for j in range( len(mjjBins)-1):
       x=mjjBins[j+1]-mjjBins[j];
       binsM.Fill(mjjBins[j]+0.5*x,x);

# masses 
h1=TH1D("Mjj_"+ka,"Mjj_"+ka, len(mjjBins)-1, mjjBins )
h2=TH1D("Mbb_"+ka,"Mbb_"+ka, len(mjjBins)-1, mjjBins)
h3=TH1D("Mjb_"+ka,"Mjb_"+ka, len(mjjBins)-1, mjjBins)
h4=TH1D("Mee_"+ka,"Mee_"+ka,200,1,401)
h5=TH1D("Mmm_"+ka,"Mmm_"+ka,200,1,401)
h6=TH1D("Mje_"+ka,"Mje_"+ka, len(mjjBins)-1, mjjBins)
h7=TH1D("Mjm_"+ka,"Mjm_"+ka, len(mjjBins)-1, mjjBins)
h8=TH1D("Mjg_"+ka,"Mjg_"+ka, len(mjjBins)-1, mjjBins)
h9=TH1D("Mgg_"+ka,"Mgg_"+ka, len(mjjBins)-1, mjjBins)
h10=TH1D("Mge_"+ka,"Mge_"+ka,len(mjjBins)-1, mjjBins)
h11=TH1D("Mbg_"+ka,"Mbg_"+ka,len(mjjBins)-1, mjjBins)
h12=TH1D("Mgm_"+ka,"Mgm_"+ka,len(mjjBins)-1, mjjBins)
h13=TH1D("Mem_"+ka,"Mem_"+ka,100,1,1001)
h14=TH1D("Mbe_"+ka,"Mbe_"+ka,len(mjjBins)-1, mjjBins)
h15=TH1D("Mbm_"+ka,"Mbm_"+ka,len(mjjBins)-1, mjjBins)


# before NN cut 
h1b=TH1D("Mjj_b_"+ka,"Mjj_b_"+ka, len(mjjBins)-1, mjjBins )
h2b=TH1D("Mbb_b_"+ka,"Mbb_b_"+ka, len(mjjBins)-1, mjjBins)
h3b=TH1D("Mjb_b_"+ka,"Mjb_b_"+ka, len(mjjBins)-1, mjjBins)
h4b=TH1D("Mee_b_"+ka,"Mee_b_"+ka,200,1,401)
h5b=TH1D("Mmm_b_"+ka,"Mmm_b_"+ka,200,1,401)
h6b=TH1D("Mje_b_"+ka,"Mje_b_"+ka, len(mjjBins)-1, mjjBins)
h7b=TH1D("Mjm_b_"+ka,"Mjm_b_"+ka, len(mjjBins)-1, mjjBins)
h8b=TH1D("Mjg_b_"+ka,"Mjg_b_"+ka, len(mjjBins)-1, mjjBins)
h9b=TH1D("Mgg_b_"+ka,"Mgg_b_"+ka, len(mjjBins)-1, mjjBins)
h10b=TH1D("Mge_b_"+ka,"Mge_b_"+ka,len(mjjBins)-1, mjjBins)
h11b=TH1D("Mbg_b_"+ka,"Mbg_b_"+ka,len(mjjBins)-1, mjjBins)
h12b=TH1D("Mgm_b_"+ka,"Mgm_b_"+ka,len(mjjBins)-1, mjjBins)
h13b=TH1D("Mem_b_"+ka,"Mem_b_"+ka,100,1,1001)
h14b=TH1D("Mbe_b_"+ka,"Mbe_b_"+ka,len(mjjBins)-1, mjjBins)
h15b=TH1D("Mbm_b_"+ka,"Mbm_b_"+ka,len(mjjBins)-1, mjjBins)

debug=TH1D("debug_"+ka,"debug_"+ka,10,0,10)

######## scaled by CM energy in GeV 
CMS_E_GEV=CMS_RUN2 # in  GeV 

# masses 
ka=ka+"_cms"

mjjBinsE=[]
for j in range( len(mjjBins)):
       mjjBinsE.append(mjjBins[j]/CMS_E_GEV) 
mjjBinsE = array("d", mjjBinsE)

binsEM = TH1D("bins_m_"+ka, "bins_m_"+ka, len(mjjBinsE)-1, mjjBinsE);
for j in range( len(mjjBinsE)-1):
       x=mjjBinsE[j+1]-mjjBinsE[j];
       binsEM.Fill(mjjBinsE[j]+0.5*x,x);


h1E=TH1D("Mjj_"+ka,"Mjj_"+ka, len(mjjBinsE)-1, mjjBinsE )
h2E=TH1D("Mbb_"+ka,"Mbb_"+ka, len(mjjBinsE)-1, mjjBinsE)
h3E=TH1D("Mjb_"+ka,"Mjb_"+ka, len(mjjBinsE)-1, mjjBinsE)
h4E=TH1D("Mee_"+ka,"Mee_"+ka,200,1,401/CMS_E_GEV)
h5E=TH1D("Mmm_"+ka,"Mmm_"+ka,200,1,401/CMS_E_GEV)
h6E=TH1D("Mje_"+ka,"Mje_"+ka, len(mjjBinsE)-1, mjjBinsE)
h7E=TH1D("Mjm_"+ka,"Mjm_"+ka, len(mjjBinsE)-1, mjjBinsE)
h8E=TH1D("Mjg_"+ka,"Mjg_"+ka, len(mjjBinsE)-1, mjjBinsE)
h9E=TH1D("Mgg_"+ka,"Mgg_"+ka, len(mjjBinsE)-1, mjjBinsE)
h10E=TH1D("Mge_"+ka,"Mge_"+ka,len(mjjBinsE)-1, mjjBinsE)
h11E=TH1D("Mbg_"+ka,"Mbg_"+ka,len(mjjBinsE)-1, mjjBinsE)
h12E=TH1D("Mgm_"+ka,"Mgm_"+ka,len(mjjBinsE)-1, mjjBinsE)
h13E=TH1D("Mem_"+ka,"Mem_"+ka,100,1,1001/CMS_E_GEV)
h14E=TH1D("Mbe_"+ka,"Mbe_"+ka,len(mjjBinsE)-1, mjjBinsE)
h15E=TH1D("Mbm_"+ka,"Mbm_"+ka,len(mjjBinsE)-1, mjjBinsE)

# before NN cut 
h1Eb=TH1D("Mjj_b_"+ka,"Mjj_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h2Eb=TH1D("Mbb_b_"+ka,"Mbb_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h3Eb=TH1D("Mjb_b_"+ka,"Mjb_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h4Eb=TH1D("Mee_b_"+ka,"Mee_b_"+ka,200,1,401/CMS_E_GEV)
h5Eb=TH1D("Mmm_b_"+ka,"Mmm_b_"+ka,200,1,401/CMS_E_GEV)
h6Eb=TH1D("Mje_b_"+ka,"Mje_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h7Eb=TH1D("Mjm_b_"+ka,"Mjm_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h8Eb=TH1D("Mjg_b_"+ka,"Mjg_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h9Eb=TH1D("Mgg_b_"+ka,"Mgg_b_"+ka, len(mjjBinsE)-1, mjjBinsE)
h10Eb=TH1D("Mge_b_"+ka,"Mge_b_"+ka,len(mjjBinsE)-1, mjjBinsE)
h11Eb=TH1D("Mbg_b_"+ka,"Mbg_b_"+ka,len(mjjBinsE)-1, mjjBinsE)
h12Eb=TH1D("Mgm_b_"+ka,"Mgm_b_"+ka,len(mjjBinsE)-1, mjjBinsE)
h13Eb=TH1D("Mem_b_"+ka,"Mem_b_"+ka,100,1,1001/CMS_E_GEV)
h14Eb=TH1D("Mbe_b_"+ka,"Mbe_b_"+ka,len(mjjBinsE)-1, mjjBinsE)
h15Eb=TH1D("Mbm_b_"+ka,"Mbm_b_"+ka,len(mjjBinsE)-1, mjjBinsE)

# write AR
anomr=TH1D("AR","Anomaly Region",4,0,4)
anomr.Fill(1,CutOutlierMC)


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



## scaled by CME
h1E.Sumw2();
h2E.Sumw2();
h3E.Sumw2();
h4E.Sumw2();
h5E.Sumw2();
h6E.Sumw2();
h7E.Sumw2();
h8E.Sumw2();
h9E.Sumw2();
h10E.Sumw2();
h11E.Sumw2();
h12E.Sumw2();
h13E.Sumw2();
h14E.Sumw2();
h15E.Sumw2();

h1Eb.Sumw2();
h2Eb.Sumw2();
h3Eb.Sumw2();
h4Eb.Sumw2();
h5Eb.Sumw2();
h6Eb.Sumw2();
h7Eb.Sumw2();
h8Eb.Sumw2();
h9Eb.Sumw2();
h10Eb.Sumw2();
h11Eb.Sumw2();
h12Eb.Sumw2();
h13Eb.Sumw2();
h14Eb.Sumw2();
h15Eb.Sumw2();

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
# scale in GeV
CMS_E=CMS_RUN2


for i in range(len(proc)):
   ev=0

   # initialize
   RMM = np.zeros(shape=(evtINchunk, mSize*mSize))
   masses=np.zeros(shape=(evtINchunk, MaxNumberOfMasses))
   weights=[]

   # total events in file
   NtotInFile=(rfile[i].inputNN).GetEntries()
   print("Analyse ",NtotInFile," from file=",rfile[i])

   for event in rfile[i].inputNN:

       Tevent=event.event

#       if (onlyFraction<1.0):
#          ran=random.uniform(0, 1)
#          if (ran>onlyFraction): continue

       if (onlyFraction<1.0):
          takeevent=False;
          eventShifted=Tevent+1 # avoid overlap with trained samples in CSV file 
          if (eventShifted%int(1/onlyFraction) == 0): takeevent=True
          if (takeevent == False): continue

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

       # scale in GeV
       CMS_E=CMS_RUN2
       # switch to a different CMS energy for Run3
       if (Trun>RUN3_MIN_RUN): CMS_E=CMS_RUN3

       # fetch RMM values
       for i3 in range(NN):
              w=inx1[i3];
              h=inx2[i3];
              val=float(a[i3])
              emptyMatrix[w][h] = val
              # mjj
              if (h==mjj[0] and w==mjj[1]): v_mjj=val 
              # mbb
              if (h==mbb[0] and w==mbb[1]): v_mbb=val 
              # mee
              if (h==mee[0] and w==mee[1]): v_mee=val 
              #mumu
              if (h==mmumu[0] and w==mmumu[1]): v_mmm=val 
              # mjb 
              if (h==mbj[0] and w==mbj[1]):  v_mjb=val 
              # mje 
              if (h==mje[0] and w==mje[1]):  v_mje=val 
              # mje 
              if (h==mjmu[0] and w==mjmu[1]):  v_mjm=val 
              # mjg 
              if (h==mjg[0] and w==mjg[1]):  v_mjg=val 
              # mbe 
              if (h==mbe[0] and w==mbe[1]):  v_mbe=val
              # mje 
              if (h==mbmu[0] and w==mbmu[1]):  v_mbm=val
              # mjg 
              if (h==mbg[0] and w==mbg[1]):  v_mbg=val


       # fill masses before any cut
       h1b.Fill( v_mjj*CMS_E, weight )
       h2b.Fill( v_mbb*CMS_E, weight )
       h3b.Fill( v_mjb*CMS_E, weight )
       h4b.Fill( v_mee*CMS_E, weight )
       h5b.Fill( v_mmm*CMS_E, weight )
       # additional
       h6b.Fill( v_mje*CMS_E, weight ) # j+e 
       h7b.Fill( v_mjm*CMS_E, weight ) # j+mu 
       h8b.Fill( v_mjg*CMS_E, weight ) # j+gamma 
       h14b.Fill( v_mbe*CMS_E, weight ) # b+e 
       h15b.Fill( v_mbm*CMS_E, weight ) # b+mu 
       h11b.Fill( v_mbg*CMS_E, weight ) # b+gamma 


       ## after division by 1/sqrt(CM) in GeV. Directly from RMM 
       h1Eb.Fill( v_mjj, weight )
       h2Eb.Fill( v_mbb, weight )
       h3Eb.Fill( v_mjb, weight )
       h4Eb.Fill( v_mee, weight )
       h5Eb.Fill( v_mmm, weight )
       # additional
       h6Eb.Fill( v_mje, weight ) # j+e 
       h7Eb.Fill( v_mjm, weight ) # j+mu 
       h8Eb.Fill( v_mjg, weight ) # j+gamma 
       h14Eb.Fill( v_mbe, weight ) # b+e 
       h15Eb.Fill( v_mbm, weight ) # b+mu 
       h11Eb.Fill( v_mbg, weight ) # b+gamma 


       # print(v_mje,  v_mjm) 

       # flatten
       dataRMM=(emptyMatrix.flatten()).tolist()
       RMM[evt,:]=dataRMM       
       masses[evt,:]=numpy.array([v_mjj,v_mbb,v_mjb,v_mee,v_mmm,v_mje,v_mjm,v_mjg,v_mbe,v_mbm,v_mbg])
       weights.append(weight) 

       evt=evt+1   # events in chunk 
       ev=ev+1     # events in this file 
       ntot=ntot+1 # all events 
       if (ev ==  NtotInFile or evt%evtINchunk==0):

                     df = pd.DataFrame(data=RMM, columns=columnsX)
                     # in the case of parameterized NN, need to add extra parameter..
                     if (NN_TYPE == "PAE"): df['Parameter'] = input_trigger_type;  # parameter that tells which trigger stream 
                     # one hot encoder
                     if (NN_TYPE.find("HAE")>-1):
                             for trg in range(1,8):  # channels 1-7 are 0 
                                  df[f'Parameter{trg}'] = 0
                             df[f'Parameter{input_trigger_type}'] = 1

                     df=df.drop(col0, axis = 1)
                     #print("Apply scalers and remove 0 columns: DF size=",df.size," DF shape=",df.shape," DF dimension=",df.ndim)
                     RMM_T = df.to_numpy()

                     if (IsStandard):
                       RMM_T = scalerStandard.transform(RMM_T)
                       RMM_T = scalerMinMax.transform(RMM_T)

                     predictions = loaded_model.predict( RMM_T )
                     if NN_TYPE == "PAE":
                       nparam = 1
                       train_loss = tf.keras.losses.mae(predictions[:, :-nparam], RMM_T[:, :-nparam]).numpy()
                     elif (NN_TYPE.find("HAE")>-1):
                       nparam = 7
                       train_loss = tf.keras.losses.mae(predictions[:, :-nparam], RMM_T[:, :-nparam]).numpy()
                     else: # standard AE 
                       train_loss = tf.keras.losses.mae(predictions, RMM_T).numpy()

                     nle=len(train_loss)
                     if (ev ==  NtotInFile):
                                          nle=ev- nchunk*evtINchunk 
                                          print(" -> Last event ",ev, " from ",NtotInFile, "remaining=",nle)

                     for ch in range(nle):

                            xloss=train_loss[ch]
                            we=weights[ch]
                            xlog= math.log(xloss) 
                            #print(xlog)
                            h_loss.Fill(xlog, we)
                            debug.Fill(1)

                            if (xlog < CutOutlierMC): continue # reject SM using MC outlier 
                            # This is side-band control region -10.10 - 10.0
                            if (workingPoint==-1):
                                               if (xlog >  CutOutlierMC+0.10): continue  

                            debug.Fill(2)
                            ma=masses[ch]
                            # after SM reject 
                            h1.Fill( ma[0]*CMS_E, we )
                            h2.Fill( ma[1]*CMS_E, we )
                            h3.Fill( ma[2]*CMS_E, we )
                            h4.Fill( ma[3]*CMS_E, we )
                            h5.Fill( ma[4]*CMS_E, we )
                            # additional masses
                            h6.Fill( ma[5]*CMS_E, we )
                            h7.Fill( ma[6]*CMS_E, we  )
                            h8.Fill( ma[7]*CMS_E, we )
                            h14.Fill( ma[8]*CMS_E, we )
                            h15.Fill( ma[9]*CMS_E, we )
                            h11.Fill( ma[10]*CMS_E, we )


                            # division by CM in GEV energy. Taken from RMM
                            # after SM reject
                            h1E.Fill( ma[0], we )
                            h2E.Fill( ma[1], we )
                            h3E.Fill( ma[2], we )
                            h4E.Fill( ma[3], we )
                            h5E.Fill( ma[4], we )
                            # additional masses
                            h6E.Fill( ma[5], we )
                            h7E.Fill( ma[6], we  )
                            h8E.Fill( ma[7], we )
                            h14E.Fill( ma[8], we )
                            h15E.Fill( ma[9], we )
                            h11E.Fill( ma[10], we )


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
  debug.Write()
  hfile.Close()


hfile=TFile(rootfile,"RECREATE","signatures")
debug.Write()
h_loss.Write()
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

## division by CMS in TeV
h1E.Write()
h2E.Write()
h3E.Write()
h4E.Write()
h5E.Write()
h6E.Write()
h7E.Write()
h8E.Write()
h9E.Write()
h10E.Write()
h11E.Write()
h12E.Write()
h13E.Write()
h14E.Write()
h15E.Write()
# before
h1Eb.Write()
h2Eb.Write()
h3Eb.Write()
h4Eb.Write()
h5Eb.Write()
h6Eb.Write()
h7Eb.Write()
h8Eb.Write()
h9Eb.Write()
h10Eb.Write()
h11Eb.Write()
h12Eb.Write()
h13Eb.Write()
h14Eb.Write()
h15Eb.Write()





# bins
binsM.Write()
binsEM.Write()
anomr.Write()

# tree
if save_tree:
    tree.Write()

hfile.Close()
print ("Write=",rootfile) 






