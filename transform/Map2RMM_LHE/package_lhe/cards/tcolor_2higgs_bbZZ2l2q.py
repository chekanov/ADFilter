# HepSim Pythia
# SM diHiggs production, decay to multi-lepton
# Similar to: https://gitlab.cern.ch/atlas-physics/pmg/mcjoboptions/-/blob/master/600xxx/600863/mc.PhPy8EG_NNPDF30NLO_HHML_bbZZ2l2q_1LcHHH01d0.py
ApplyParticleSlim=off
# do not apply slim
#
#numbers of events 
EventsNumber=100000
#Collision settings
Random:setSeed = on
Random:seed = 0
Beams:idA = 2212
Beams:idB = 2212
Beams:eCM = 13000.
#physics processes
Beams:frameType = 4
HardQCD:all = off
#
ParticleDecays:tau0Max = 10
Tune:pp = 14
Tune:ee = 7

# PDF:useLHAPDF = on
SpaceShower:rapidityOrder = on
SigmaProcess:alphaSvalue = 0.140
SpaceShower:pT0Ref = 1.56
SpaceShower:pTmaxFudge = 0.91
SpaceShower:pTdampFudge = 1.05
SpaceShower:alphaSvalue = 0.127
TimeShower:alphaSvalue = 0.127
BeamRemnants:primordialKThard = 1.88
MultipartonInteractions:pT0Ref = 2.09
MultipartonInteractions:alphaSvalue = 0.126

#Makes particles with c*tau>10 mm stable
ParticleDecays:tau0Max = 10
# Z decays
25:oneChannel = on 0.5   100 23 23
# b-decays
25:addChannel = on 0.5   100 5 -5   # bb decay 
23:onMode = off
23:onIfAny= 1 2 3 4 5 11 13 15
# W/Z minimum mass
24:mMin = 0
24:mMax = 99999
23:mMin = 0
23:mMax = 99999
TimeShower:mMaxGamma = 0
