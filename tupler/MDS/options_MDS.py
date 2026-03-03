from PyConf.reading import get_particles, get_pvs, get_mc_particles
from DaVinci import make_config, Options
from DaVinci.algorithms import create_lines_filter
from FunTuple import FunTuple_Particles as Funtuple
import Functors as F
from Hlt2Conf.lines.bandq.builders import c_hadrons, charged_hadrons, b_hadrons
from Hlt2Conf.lines.charmonium_to_dimuon import make_jpsi
from GaudiKernel.SystemOfUnits import MeV, picosecond, mm
from RecoConf.algorithms_thor import ParticleContainersMerger, ParticleCombiner
from FunTuple import FunTuple_MCParticles as FuntupleMC
from RecoConf.event_filters import require_pvs

from RecoConf.muonid import make_muon_hits
from PyConf.Algorithms import PrintDecayTree
from RecoConf import standard_particles

from PyConf.Algorithms import PrintMCTree

import utils


def diyline():
    
    muons = standard_particles.make_long_muons()
    #muons = standard_particles.make_long_cb_muons()
    #muons = standard_particles.make_down_muons()
    

    combination_code = F.require_all(
        F.MASS > -0.1,
    )

    vertex_code = F.require_all(
        F.MASS > -0.1,
    )

    dimuon = ParticleCombiner(
        name="MyDiMuon",
        Inputs=[muons, muons],
        DecayDescriptor='J/psi(1S) -> mu+ mu-',
        CombinationCut=combination_code,
        CompositeCut=vertex_code)
    return dimuon

def make_algs_MDS(options, line, algname):
    
    data = diyline()
    line_type_name = "Turbo"
    if options.simulation:
        line_type_name = "HLT2"
    else:
        if options.input_process == "TurboSpruce":
            line_type_name = "Turbo"
        if options.input_process == "Spruce":
            line_type_name = "Spruce"

    data_path = f"/Event/{line_type_name}/{line}/Particles"        
    print("====INFO: using data path: "+data_path)
    data = get_particles(data_path)
    line_prefilter = create_lines_filter(name=f"PreFilter_{line}", lines=[line])

    pdt = PrintDecayTree(name="PrintLine"+line, Input=data)

    pvs = get_pvs()

    fields = {
        "Jpsi" :   "J/psi(1S) -> mu+ mu-",
        "Mup"  :   "J/psi(1S) -> ^mu+ mu-",
        "Mum"  :   "J/psi(1S) -> mu+ ^mu-",
    }


    all_vars = utils.make_base_vars(pvs)
    #all_vars += utils.make_Hlt2TisTos_vars(data)
    if options.simulation == True:
        all_vars += utils.make_MCTruth_vars(data, algname)
    vertex_vars = utils.make_vertex_vars(pvs)
    vertex_vars_3daughters = utils.make_vertex_vars_3daughters()
    vertex_vars_4daughters = utils.make_vertex_vars_4daughters()
    track_vars = utils.make_track_vars(pvs)

    variables = {
        "ALL": all_vars,
        #"Jpsi": vertex_vars,
        #"Mum": track_vars,
        #"Mup": track_vars,
    }

    evt_vars = utils.make_evt_vars(pvs, line, "Hlt2", options.simulation)

    funtuple = Funtuple(
        name=algname,
        tuple_name="DecayTree",
        fields=fields,
        variables=variables,
        event_variables=evt_vars,
        inputs=data,
    )

    mct_fields = {
        "Bcp": "[B_c+ -> (H_10 -> pi+ mu-) mu+]CC",
        "HNL": "[B_c+ -> ^(H_10 -> pi+ mu-) mu+]CC",
        "PipHNL": "[B_c+ -> (H_10 -> ^pi+ mu-) mu+]CC",
        "MumHNL": "[B_c+ -> (H_10 -> pi+ ^mu-) mu+]CC",
        "Mup": "[B_c+ -> (H_10 -> pi+ mu-) ^mu+]CC",
    }

    mc_particles = get_mc_particles('/Event/MC/Particles')

    printMC = PrintMCTree(
        MCParticles=mc_particles, ParticleNames=["B_c+", "B_c-"], OutputLevel=4
    )


    variables_mc = utils.make_mc_variables()
    mct_variables = { 'ALL': variables_mc }
    funtuple_mc = FuntupleMC(
        name=algname+"MC",
        tuple_name="DecayTree",
        fields=mct_fields,
        variables=mct_variables,
        inputs=mc_particles,
    )

    if options.simulation == True:
        #algs = [line_prefilter, require_pvs(pvs), funtuple]
        #algs = [require_pvs(pvs), funtuple]
        algs = [funtuple]
        #algs = [printMC]
        #algs = [line_prefilter, pdt]
        return algs
    else:
        #algs = [line_prefilter, require_pvs(pvs), funtuple]
        #algs = [require_pvs(pvs), funtuple]
        #algs = [funtuple]
        algs = [line_prefilter, pdt]
        return algs

def main(options: Options):
    #line = "Hlt2QEE_MDS_BDT_nHits"
    #line = "Hlt2QEE_JpsiToMuMu_Detached"
    #line = "Hlt2QEE_SingleHighPtMuonFull"
    #line = "Hlt2QEE_SingleHighPtMuonIsoFull"
    line = "Hlt2BandQ_DiMuonIncFull"
    #line = "Hlt2BandQ_DiMuonSoftFull"
    #line = "Hlt2QEE_DisplacedDiMuon"

    #line = "SpruceQEE_ZToMuMu"
    #line = "SpruceQEE_SingleJet15"
    #line = "Hlt2QEE_JpsiToMuMu_Prompt"
    #line = "SpruceQEE_ZToMuMu"
    name = "myTuple"
    algs = {name : make_algs_MDS(options, line, name) }

    return make_config(options, algs)
