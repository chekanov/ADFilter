import random
import sys
sys.path.append("modules/")

from array import *
from math import *
import math,sys,os 
from array import array
from decimal import Decimal
import numpy
import random
import sys,zipfile,json,math

from array import *
from math import *

CMS=13000.0

# Run 2 and Run3 in TeV
CMS_RUN2=CMS
CMS_RUN3=13600.0

# Run2 and Run3 in TeV
CMS_RUN2_TEV=CMS_RUN2*0.001
CMS_RUN3_TEV=CMS_RUN3*0.001

# minimum Run for Run2. After this run number, run3 starts
RUN3_MIN_RUN=420000


lumi2015=3244.54  # error 3.21 +-0.07 
Lumi2015=" %.1f fb^{-1}" % (lumi2015/1000.)
intLUMI2015="#int L dt = "+Lumi2015

lumi2016=33402.2  # error 3.21 +-0.07 
Lumi2016=" %.1f fb^{-1}" % (lumi2016/1000.)
intLUMI2016="#int L dt = "+Lumi2016

lumi2017=44630.6   # error 3.21 +-0.07 
Lumi2017=" %.1f fb^{-1}" % (lumi2017/1000.)
intLUMI2017="#int L dt = "+Lumi2017

# Nov
# lumi2018=43003.7 # error 3.21 +-0.07 
# Dec
lumi2018=58791.6 # error 3.21 +-0.07 
Lumi2018=" %.1f fb^{-1}" % (lumi2018/1000.)
intLUMI2018="#int L dt = "+Lumi2018

# lumi  in pb
lumi2015_2018=lumi2015+lumi2016+lumi2017+lumi2018
lumi=lumi2015_2018 # take into account missing files 

# new recomendation
# https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LuminosityForPhysics#2015_2018_13_TeV_proton_proton_f

Lumi=" %.0f fb^{-1}" % (lumi/1000.)
intLUMI="#int L dt = "+Lumi

Lumi10=" %.0f fb^{-1}" % (0.1*lumi/1000.)
intLUMI10="#int L dt = "+Lumi10


#########################################################
# cut to select outlier events

# 20 pb WP
CutOutlier_20PB=-9.39

# MC region for 10 pb working point ("data limit") 
CutOutlier_10PB=-9.10 

# WP region for 1 pb working point 
CutOutlier_1PB=-8.0

# WP region 0.1 pb
CutOutlier_01PB=-6.50 


# main cuts
CutOutlierData=CutOutlier_01PB
CutOutlierMC=CutOutlier_10PB

# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/XsecSummaryWjetsPowPy8Incl
# 361100-361104
# 2.65 + 2.65  + 4.83 + 2.41
powheg_lumi=(2.65+2.65+4.83+2.41)*1000  # pb-1
powheg_kfactor=1.0172
powheg_scale=(lumi/powheg_lumi)*powheg_kfactor

# https://twiki.cern.ch/twiki/bin/view/AtlasProtected/XsecSummaryTTbar
ttbar_cross=729.0 # pb
# 97% grid efficiency
ttbar_events=49874000*0.97
# ttbar_lumi=ttbar_file.Get("cutflow").GetBinContent(1) /(ttbar_cross*0.543) # pb
ttbar_lumi=ttbar_events /(ttbar_cross*0.543) # pb
ttbar_kfactor=1.195
ttbar_scale=(lumi/ttbar_lumi)*ttbar_kfactor

# this is ROOT file after using the model with 1% (nominal)
root_model2use="root/analysis_median.root"

# https://twiki.cern.ch/twiki/bin/view/AtlasProtected/XsecSummarySingleTop
# Run 410659 and  410647 
stop_cross=22.175+36.996 # pb
# stop_events=5968000+6226000+6226000+6226000.. Rough number of input 
stop_events=9968000*6
# ttbar_lumi=ttbar_file.Get("cutflow").GetBinContent(1) /(ttbar_cross*0.543) # pb
stop_lumi=stop_events /stop_cross # pb
stop_kfactor=1.195
stop_scale=(lumi/stop_lumi)*stop_kfactor


# Dijet QCD
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissMC15
# in nb-1
pythia_lumi={}

# events/cros*eff 
# in nb
factors={}
factors[0]=2000000 / (7.8420E+07*1.0240E+0)
factors[1]=2000000 / (7.8420E+07*6.7198E-04)
factors[2]=1992000  / (2.4334E+06*3.3264E-04)
factors[3]=1767000  / (2.6454E+04*3.1953E-04)
factors[4]=1997000 /  (2.5464E+02*5.3009E-04)
factors[5]=1995000  / (4.5536E+00*9.2325E-04)
factors[6]=1997000 / (2.5752E-01*9.4016E-04)
factors[7]=1990000  / (1.6214E-02*3.9282E-04)
factors[8]=2000000 / (6.2505E-04*1.0162E-02)
factors[9]=2000000 / (1.9640E-05*1.2054E-02)
factors[10]=2000000 / (1.1961E-06*5.8935E-03)



DataLab="Data #sqrt{s}=13 TeV"
KinemCuts="E_{T}^{jet}>410 GeV  |#eta^{#gamma}|<2.5";
ATLASprel="ATLAS internal"
mcBkg="PYTHIA8"
mcSig="PYTHIA t#bar{t}"
mcPowheg="W+jet POWHEG"
mcPowhegTTbar="t#bar{t} POWHEG"
mcSTop="s-top POWHEG"
mcBkgHrw="HERWIG++ QCD"
mcBkgWJ="Multijets PYTHIA8"

DataLab2015="Data 2015 #sqrt{s}=13 TeV"
DataLab2016="Data 2016 #sqrt{s}=13 TeV"
DataLab2017="Data 2017 #sqrt{s}=13 TeV"
DataLab2018="Data 2018 #sqrt{s}=13 TeV"


mjjBinsL = [99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156];

mjjBins = array("d", mjjBinsL)


from os.path import exists

confile='data/config.json'
file_exists = exists( confile )

if (file_exists):
 with open( confile ) as json_file:
    data = json.load(json_file)
    maxNumber=int(data['maxNumber'])
    maxTypes=int(data['maxTypes'])
    mSize=int(data['mSize'])
    print("Read from file ",confile)
else:
    maxNumber=10
    maxTypes=5
    mSize=5


print ("maxNumber=",maxNumber," maxTypes=",maxTypes," mSize=",mSize) 
mSize=maxTypes*maxNumber+1;

######################### define position of invarinat masses from RMM #################

# dijet invariant mass
x=1+0*maxNumber+1  # X position  
y=1+0*maxNumber    # Y position 
mjj=(x,y) #  index of Mjj  matrix ellement 

# PT of first jet
x=1+0*maxNumber  # X position  
y=1+0*maxNumber  # Y position 
pT=(x,y) #  index of Mjj  matrix ellement 

#  bb mass 
x=1+1*maxNumber+1
y=1+1*maxNumber
mbb=(x,y)

#  bj mass 
x=1+1*maxNumber
y=1+0*maxNumber
mbj=(x,y) 

# mu+mu 
x=1+2*maxNumber+1
y=1+2*maxNumber
mmumu=(x,y)

# e+e 
x=1+3*maxNumber+1
y=1+3*maxNumber
mee=(x,y)

# j+mu 
x=1+2*maxNumber
y=1+0*maxNumber
mjmu=(x,y)

# j+e 
x=1+3*maxNumber
y=1+0*maxNumber
mje=(x,y)

# j+gamma 
x=1+4*maxNumber
y=1+0*maxNumber
mjg=(x,y)

# b+mu 
x=1+2*maxNumber
y=1+1*maxNumber
mbmu=(x,y)

# b+e  
x=1+3*maxNumber
y=1+1*maxNumber
mbe=(x,y)

# b+gamma 
x=1+4*maxNumber
y=1+1*maxNumber
mbg=(x,y)

############# end invariant mass definitions using RMM ############

### This list contains excluded values for Z-score calculation
### We excluding pT of leading jet, Mjj and mbb
# excluded_val= ( pT, mjj, mbb)
# excluded_val= (mjj, mbb)
# print ("Excluded cells=",excluded_val ) 

#### Exclusion values for RMM matrix #############
###################################


# dijet invariant mass
x=2 # X position  
y=1 # Y position 
inx1=x*mSize+y; #  index of hxw matrix ellement 

# pT1 
x=1 # X position  
y=1 # Y position 
inx2=x*mSize+y; #  index of hxw matrix ellement 

# Mjj for for light-jet + b-jets
x=1+maxNumber # X position  
y=1 # Y position 
inx3=x*mSize+y; #  index of hxw matrix ellement 

# pT1 for for b-jets
x=1+maxNumber # X position  
y=1+maxNumber # Y position 
inx4=x*mSize+y; #  index # pT for for b-jets

# Mjj for for 2-b jets 
x=2+maxNumber # X position  
y=1+maxNumber # Y position 
inx5=x*mSize+y; #  index of hxw matrix ellement 

mjjBinsL = [99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156];

mjjBins = array("d", mjjBinsL)



# exlusion matrix for RMM in terms of indexes (how it is packed) 
excluded=(inx1,inx2,inx3,inx4,inx5)

# Save mathplot in CSV file
import csv
def SavePlotXY(xfile,lines, Xlab="X", Ylab="Y"):
   #NrLines=len(lines)
   #print("Nr of lines",NrLines)
   print("Save plot in CSV ",xfile);
   with open(xfile, 'w') as myfile:
            data=lines[0].get_data()
            writer = csv.writer(myfile)
            writer.writerow([Xlab, Ylab])
            for i in range(len(data[0])):
                writer.writerow([data[0][i], data[1][i]])

from numpy import savetxt
def SaveNumpyData(xfile,lines):
    savetxt(xfile,lines,delimiter=',')


# SSM
WPrime=[
         [801123,750,4.50E-03*1000,10000,4289,3756],
         [801124,1250,3.41E-04*1000,10000,4833,4397],
         [801125,2250,1.07E-05*1000,10000,5471,5111],
         [801126,3250,7.27E-07*1000,10000,5766,5437],
         [801127,4250,6.49E-08*1000,10000,6090,5732],
         [801128,5250,5.51E-09*1000,10000,6205,5833],
         [801129,6250,4.61E-10*1000,10000,6359,6024]
       ]

# Nr of events for this lumi
WPrimeEvents={}
for j in range(len( WPrime )):
       mmm=WPrime[j]
       ngen=float(mmm[3])
       effacc=(mmm[4]/ngen) * (mmm[5]/ngen) # eff*accep
       effacc=1
       WPrimeEvents[mmm[0]] = mmm[2]*lumi * effacc
  

########### Simplified DM models
ffc=1 # convert nb to pb (already in pb?!) 
# run,mass,xcross(pb),events, Mjjll, Mjjl accepted 
WZsim=[
         [502946,700,0.1119509*ffc,20000, 824,  9799],
         [502947,1200,0.024642*ffc,20000, 1300, 11997],
         [502948,2000,0.003264*ffc,20000, 1688, 13402],
         [502949,3000,0.0003579*ffc,20000, 1866, 14050],
         [502950,4000,0.0000463*ffc,20000, 1988, 14521],
         [502951,5000,0.000006583*ffc,20000, 2147, 14804],
         [502952,6000,0.000001031*ffc,20000, 2285, 14827]
        ]

# Nr of events for this lumi
DMsimEvents={}
for j in range(len( WZsim )):
       mmm=WZsim[j]
       ngen=float(mmm[3])
       effacc=(mmm[4]/ngen) * (mmm[5]/ngen) # eff*accep
       effacc=1
       DMsimEvents[mmm[0]] = mmm[2]*lumi * effacc  # events*eff*acc 


###### Radion model #################### 
ffc=1 # convert nb to pb (already in pb?!) 
# run,mass,xcross(pb),events, Mjjll, Mjjl accepted 
RadionM=[
         [506939,500, 2.26*ffc,20000, 14,  7961],
         [506940,1000, 0.0658*ffc,20000, 145, 10438],
         [506941,2000, 0.0005266*ffc,20000, 463, 11769],
         [506942,3000, 0.0000153837*ffc,20000, 813, 12493],
         [506943,4000, 0.0000006825*ffc,20000, 1147, 13106],
         [506944,5000, 0.0000000305549*ffc,20000, 1430, 13850],
         [506945,6000, 0.0000000012557*ffc,20000, 1767, 14454]
        ]

# Nr of events for this lumi
DRadionEvents={}
for j in range(len( RadionM  )):
       mmm=RadionM[j]
       ngen=float(mmm[3])
       effacc=(mmm[4]/ngen) * (mmm[5]/ngen) # eff*accep
       effacc=1
       DRadionEvents[mmm[0]] = mmm[2]*lumi * effacc  # events*eff*acc 

############ Composite Lepton ##########################################
ffc=1 # convert nb to pb (already in pb?!) 
# run,mass,massX,xcross(pb),events, Mjjll, Mjjl accepted 
ComLep=[
   [506946, 250, 500, 16270, 14622, 1.39E+00],
   [506947, 250, 1000, 17670, 16375, 6.05E-02],
   [506948, 500, 1000, 18254, 17729, 6.00E-02],
   [506949, 750, 1000, 18522, 17669, 2.81E-02],
   [506950, 250, 2000, 17544, 16401, 7.84E-04],
   [506951, 500, 2000, 18331, 17999, 1.07E-03],
   [506952, 750, 2000, 18484, 18231, 9.22E-04],
   [506953, 1000, 2000, 18623, 18358, 7.25E-04],
   [506954, 1250, 2000, 18488, 18215, 5.03E-04],
   [506955, 1500, 2000, 18200, 17925, 2.68E-04],
   [506956, 1750, 2000, 17850, 17343, 8.56E-05],
   [506957, 500, 3000, 18187, 17680, 5.83E-05],
   [506958, 1000, 3000, 18443, 18220, 3.40E-05],
   [506959, 1500, 3000, 18370, 18212, 2.01E-05],
   [506960, 2000, 3000, 17924, 17761, 9.05E-06],
   [506961, 2500, 3000, 17411, 17161, 2.48E-06],
   [506962, 500, 4000, 18123, 17423, 9.86E-06],
   [506963, 1000, 4000, 18455, 18140, 3.41E-06],
   [506964, 1500, 4000, 18226, 18068, 1.53E-06],
   [506965, 2000, 4000, 17947, 17819, 7.91E-07],
   [506966, 2500, 4000, 17399, 17286, 3.79E-07],
   [506967, 3000, 4000, 17060, 16893, 1.60E-07],
   [506968, 3500, 4000, 16611, 16379, 4.68E-08]
  ]


# Nr of events for this lumi
DComLep={}
for j in range(len( ComLep  )):
       mmm=ComLep[j]
       ngen=float(mmm[2])
       effacc=(mmm[3]/ngen) * (mmm[3]/ngen) # eff*accep
       effacc=1
       DComLep[mmm[0]] = mmm[5]*lumi * effacc  # events*eff*acc 



ffc=1 # convert nb to pb (already in pb!) 
# Charged Higgs H+t cross sections for tan(beta)=1
# run,mass,xcross(pb),events


# NEW 2019!
# include more mass points for mhmodm. pb for xSec*BR
# run,mass,xcross(pb),events, accepted (truth)
HPlusTtanb05_mhmodm=[
     [450597,600, 1.33184498853*ffc,56076,0,0],
     [450001,700, 0.733653118278*ffc,67642,0,0],
     [450002,800, 0.422458344*ffc,14891,0,0],
     [450003,900, 0.2509698854*ffc,13416,0,0],
     [450004,1000, 0.15352390059*ffc,18490,0,0],
     [450598,1200, 0.0615319788762*ffc,81794,0,0],
     [450599,1400, 0.0265022945008*ffc,20921,0,0],
     [450600,1600, 0.012104770414*ffc,30215,0,0],
     [450601,1800, 0.0057522212784*ffc,20921,0,0],
     [450602,2000, 0.00305527*ffc,55382,0,0]
]

DHPlusEvents={}
HPlus=[]
# recalculate input events
for j in range(len(HPlusTtanb05_mhmodm)):
       mmm=HPlusTtanb05_mhmodm[j]
       run=mmm[0]
       mass=mmm[1]
       xsec=mmm[2]
       #tt=TFile("../analysis/out/sys0/mc_hplus2021/"+str(run)+".root")
       inputevents=mmm[3]
       #inputevents=tt.Get("cutflow").GetBinContent(1)
       #print run, mass, xsec, inputevents
       effacc=1
       HPlus.append([run,mass,xsec,inputevents])
       DHPlusEvents[mmm[0]] = mmm[2]*lumi * effacc  # events*eff*acc 


# Nr of events for this lumi
WPrimeEvents={}
for j in range(len( WPrime )):
       mmm=WPrime[j]
       ngen=float(mmm[3])
       WPrimeEvents[mmm[0]] = mmm[2]*lumi *(mmm[4]/ngen) * (mmm[5]/ngen) # events*eff*acc 


def SavePlotHisto(xfile,ax):
   print("Save histogram in CSV ",xfile);
   p = ax.patches
   with open(xfile, 'w') as myfile:
            writer = csv.writer(myfile)
            writer.writerow(["Xlow", "Height"])
            for i in range(len(p) ):
                lower_left_corner=p[i].get_xy() 
                #writer.writerow([ lower_left_corner[0], p[i].get_width(), p[i].get_height()  ])
                #writer.writerow([ lower_left_corner[0], p[i].get_height()  ])
                writer.writerow([ lower_left_corner[0], lower_left_corner[1]  ])


