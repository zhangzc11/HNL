# file /home/zhicaiz/GaussDev_v56r11/Gen/DecFiles/options/42912012.py generated: Thu, 12 Feb 2026 13:27:43
#
# Event Type: 42912012
#
# ASCII decay Descriptor: pp -> ( WR- -> ( mu- (nu_Rmu -> mu- jet ) ) )
#
genAlgName="Generation"
from Gaudi.Configuration import *
importOptions( "/home/zhicaiz/GaussDev_v56r11/Gen/DecFiles/options/WRmuHeavyN.py" )
from Configurables import Generation
Generation(genAlgName).EventType = 42912012
Generation(genAlgName).SampleGenerationTool = "Special"
from Configurables import Special
Generation(genAlgName).addTool( Special )
Generation(genAlgName).Special.ProductionTool = "Pythia8Production"
from Configurables import ToolSvc
from Configurables import EvtGenDecay
ToolSvc().addTool( EvtGenDecay )
ToolSvc().EvtGenDecay.UserDecayFile = "/home/zhicaiz/GaussDev_v56r11/Gen/DecFiles/dkfiles/W_mumuqq,mN=10GeV,tN=0ps.dec"
Generation(genAlgName).Special.CutTool = ""
Generation(genAlgName).FullGenEventCutTool = "LoKi::FullGenEventCut/HeavyNFromRWInAcceptance"
from Configurables import LHCb__ParticlePropertySvc
LHCb__ParticlePropertySvc().Particles = [ "nu_Rmu     943     9900014    0       10.0    1.0e-24   unknown  9900014   0.00000" , "W_R-       947     -9900024    -1    80.38500  1.0e-24   unknown  -9900024   0.00000", "W_R+       946     9900024    1    80.38500  1.0e-24   unknown  9900024   0.00000" ]

Generation(genAlgName).Special.Pythia8Production.Commands += [
             "9900014:oneChannel = 1 0.2705271 100 13 -1 2",
             "9900014:addChannel = 1 0.0116439 100 13 -1 4",
             "9900014:addChannel = 1 0.0144031 100 13 -3 2",
             "9900014:addChannel = 1 0.2034032 100 13 -3 4",
             "9900014:addChannel = 1 0.0000006 100 13 -5 2",
             "9900014:addChannel = 1 0.0000221 100 13 -5 4",
             "9900014:addChannel = 1 0.2705271 100 -13 1 -2",
             "9900014:addChannel = 1 0.0116439 100 -13 1 -4",
             "9900014:addChannel = 1 0.0144031 100 -13 3 -2",
             "9900014:addChannel = 1 0.2034032 100 -13 3 -4",
             "9900014:addChannel = 1 0.0000006 100 -13 5 -2",
             "9900014:addChannel = 1 0.0000221 100 -13 5 -4",
             "9900024:oneChannel = 1 1 0 -13 9900014"
      ]
from Configurables import LoKi__FullGenEventCut
Generation(genAlgName).addTool( LoKi__FullGenEventCut, "HeavyNFromRWInAcceptance" )
tracksInAcc = Generation(genAlgName).HeavyNFromRWInAcceptance
tracksInAcc.Code = " count ( isGoodVfromW ) > 0 "
tracksInAcc.Preambulo += [
"from GaudiKernel.SystemOfUnits import ns, GeV, mrad"
, "isHeavyN      = ( (GDECTREE('[ X -> mu+ ... ]CC')) & ( GTHETA < 400.0*mrad ) )"
, "isRW          = ( ('W_R-' == GABSID) )"
, "isGoodMuon    = ( ( GABSID == 13 ) & (~GVEV) & ( GP > 2.0*GeV ) & ( GTHETA < 400.0*mrad ) & ( GETA > 2) )"
, "isGoodVfromN  = ( isHeavyN & ( GNINTREE( isGoodMuon, HepMC.descendants ) >= 1 ) )"
, "isGoodVfromW  = ( isRW & ( GNINTREE( isGoodVfromN, HepMC.descendants ) == 1 ) & ( GNINTREE( isGoodMuon, HepMC.descendants ) >= 1 ) )"
]

