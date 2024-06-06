# HepSim Pythia seetting
# Wprime->WZprime->lnuqq 400GeV with NNPDF23LO PDF and A14 tune
# pp>Wprime>WZprime>lnuqq
# S.Chekanoc, A.Milic 
# apply particle slim?
ApplyParticleSlim=off
# do not apply slim
# ApplyParticleSlim=off
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

