# file /home/zhicaiz/GaussDev_v56r11/Gen/DecFiles/options/42912023.py generated: Thu, 12 Feb 2026 13:27:43
#
# Event Type: 42912023
#
# ASCII decay Descriptor: pp -> (Z-> mu mu) (j j~)
#
genAlgName="Generation"
from Configurables import Generation
Generation(genAlgName).EventType = 42912023
Generation(genAlgName).SampleGenerationTool = "Special"
from Configurables import Gauss
sampleGenToolsOpts = { "Generator" : "Madgraph" }
Gauss().SampleGenerationToolOptions.update(sampleGenToolsOpts)
from Configurables import ToolSvc
from Configurables import EvtGenDecay
ToolSvc().addTool( EvtGenDecay )
ToolSvc().EvtGenDecay.UserDecayFile = "/home/zhicaiz/GaussDev_v56r11/Gen/DecFiles/dkfiles/Z_mumujj_j=udsgc_Madgraph.dec"
from Configurables import Special
Generation(genAlgName).addTool( Special )
Generation(genAlgName).Special.CutTool = ""
Generation(genAlgName).FullGenEventCutTool = "LoKi::FullGenEventCut/TwoLightTwoLeptonFromZ0"

# Configure the event type.
from Configurables import (Generation, Special)
from GaudiKernel import SystemOfUnits
from Gaudi.Configuration import importOptions

# Generation options.
Generation(genAlgName).PileUpTool           = "FixedLuminosityForRareProcess"
Generation(genAlgName).DecayTool            = ""
Generation(genAlgName).SampleGenerationTool = "Special"

# Special options.
Generation(genAlgName).addTool(Special)
Generation(genAlgName).Special.CutTool        = ""
Generation(genAlgName).Special.DecayTool      = ""

# Madgraph options.
from Configurables import Gauss
from GaudiKernel import SystemOfUnits

sampleGenToolsOpts = {
    "Commands": [ "define lj g u d s c u~ d~ s~ c~",
                  " generate p p > mu+ mu- lj lj [QCD]",
                  " set mll_sf 40" # Min invariant mass of l+l- (same flavour) lepton pair
             	  ],
    "DecEff": 0.35842 # The decfile level efficiency.
}
Gauss().SampleGenerationToolOptions.update(sampleGenToolsOpts)

# Generation cut
from Configurables import LoKi__FullGenEventCut
Generation(genAlgName).addTool( LoKi__FullGenEventCut, "TwoLightTwoLeptonFromZ0" )
tracksInAcc = Generation(genAlgName).TwoLightTwoLeptonFromZ0
tracksInAcc.Code = " ( (count ( isGoodLight ) > 1) & (count ( isGoodLepton ) > 0)) "
tracksInAcc.Preambulo += [
    "from GaudiKernel.SystemOfUnits import  GeV, mrad"
    , "isGoodLight   = ( GINTREE( (('u' == GABSID) | ('d' == GABSID) | ('s' == GABSID) | ('c' == GABSID) | ('g' == GABSID)) & ( GTHETA < 400.0*mrad )))"
    , "isGoodLepton   = ((  'Z0' == GABSID ) & GINTREE( GLEPTON & ( GTHETA < 400.0*mrad )))"
]


