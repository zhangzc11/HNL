from Gauss.Configuration import GenInit

GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 1
GaussGen.RunNumber = 1082

from Configurables import LHCbApp
LHCbApp().DDDBtag = '2024-v00.04'
LHCbApp().CondDBtag = 'sim10-2024.W31-v00.00-mu100'
LHCbApp().EvtMax = 100
