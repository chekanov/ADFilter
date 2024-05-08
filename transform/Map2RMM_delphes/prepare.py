# unpack ROOT file with RMM matrices
# you need pyROOT and numpy
# S.Chekanov (ANL)

# input file
proc=["output.root"]

from ROOT import gROOT,gPad,gStyle,TCanvas,TSpline3,TFile,TLine,TLatex,TAxis,TLegend,TPostScript
from ROOT import TH2D,TF2,TArrow,TCut,TPad,TH1D,TF1,TObject,TPaveText,TGraph,TGraphErrors,TGraphAsymmErrors
from ROOT import TGraph2D,TTree,TMultiGraph,TBranch,gSystem,gDirectory
from ROOT import TPaveStats,TProfile2D
import math
import numpy
import ROOT
import sys,math

gROOT.Reset()

name="projection"
nameX=""
nameY=""
Ymin=0.0
Ymax=500000
Xmin=0
Xmax=6.0 
ZSmin=0.00000001
ZSmax=0.5

gStyle.SetLabelSize(0.035,"xyz");

for i in range(len(proc)):
   print i, "Input file=",proc[i]
 
rfile=[]

print "Look at ", len(proc) ," processes using RMM", proc 

 
ROOT.gROOT.ProcessLine('.L Loader.C+')
for i in proc:
     rfile.append(ROOT.TFile.Open(i))
     print i

test=[]
pat=[]
validation=[]
ev=0

out=[0,0,0,0,0,0,0]

# get dimensions of the stored matrix
dimensions=(rfile[0]).Get("dimensions");
maxNumber=int(dimensions.GetBinContent(2))
maxTypes=int(dimensions.GetBinContent(3))
mSize=int(dimensions.GetBinContent(4))
print "maxNumber=",maxNumber," maxTypes=",maxTypes," mSize=",mSize
mSize=maxTypes*maxNumber+1;
print "maxNumber=",maxNumber," maxTypes=",maxTypes," mSize=",mSize

# exclude some values if needed
# dijet invariant mass
x=2 # X position  
y=1 # Y position 
inx1=x*mSize+y; #  index of hxw matrix ellement 

# pT of the leading jet 
x=1 # X position  
y=1 # Y position 
inx2=x*mSize+y; #  index of hxw matrix ellement 

# exlusion matrix
excluded=(inx1,inx2)
print "Excluded cells=",excluded

hhD = TProfile2D("profile", "profile", mSize, 0, mSize, mSize, 0, mSize, 0, 1000);
names=["MET","j", "#mu", "e", "#gamma"]
Names1=[]
Names1.append(names[0]);
for h in range(1,maxTypes+1,1):
       for i in range(1,maxNumber+1):
                 Names1.append(names[h]+"_{"+str(i)+"}");
Names2=[]
for i in range(len(Names1)):
         Names2.append(Names1[i]);
Names1= Names1[::-1]
print Names1
for h in range(mSize):
      for w in range(mSize):
        i1=h;
        i2=w;
        hhD.Fill(Names2[i1],  Names1[i2],0);
hhD.SetTitle("")
hhD.SetStats(0)
hhD.GetZaxis().SetRangeUser(ZSmin,ZSmax);
gStyle.SetPaintTextFormat(".0e");


for event in rfile[IND].inputNN:
      NN=(event.proj).size()
      a=event.proj
      inx1=event.proj_index1
      inx2=event.proj_index2
      data=[]

      emptyMatrix = numpy.zeros(shape=(mSize,mSize))
      for i3 in range(NN):
                d1=inx1[i3]; d2=inx2[i3];
                emptyMatrix[d1][d2] = float(a[i3])    # remove energy related mass cell

      # flattern
      data=(emptyMatrix.flatten()).tolist()

      #exclude some values
      for m in excluded: data[m]=0

      # validate
      for k in range(len( data )):
                          x=k/mSize
                          y=k%mSize
                          i1=x;
                          i2=mSize-y-1;
                          hhD.Fill(Names2[i1],  Names1[i2],float(data[k]));


