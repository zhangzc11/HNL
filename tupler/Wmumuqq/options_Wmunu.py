from PyConf.reading import get_particles, get_pvs, get_mc_particles
from DaVinci import make_config, Options
from DaVinci.algorithms import create_lines_filter
from FunTuple import FunTuple_Particles as Funtuple
import Functors as F
import Functors.math as fmath
from FunTuple import FunctorCollection
#from Hlt2Conf.standard_jets import make_jets, make_dijets, make_Trijets

from GaudiKernel.SystemOfUnits import MeV, picosecond, mm, GeV
from RecoConf.algorithms_thor import ParticleContainersMerger, ParticleCombiner, ParticleFilter
from FunTuple import FunTuple_MCParticles as FuntupleMC
from RecoConf.event_filters import require_pvs

from PyConf.Algorithms import PrintDecayTree
from RecoConf import standard_particles

from PyConf.Algorithms import PrintMCTree

import utils



def diyline():

    muons = ParticleFilter(
        standard_particles.make_long_muons(),
        F.FILTER(F.require_all(F.PT > 15 * GeV, F.ISMUON == True, F.INECAL == True, F.INHCAL == True, F.OWNPVIP < 0.04))
        )

    return muons
  
def make_algs_Wmunu(options, algname):
    data = diyline()

    pvs = get_pvs()

    fields = {
        "Mu" :   "[mu-]CC",
    }

    all_vars = utils.make_base_vars(pvs)
    all_vars += utils.make_HltTisTos_vars(data, Hlt1_decisions_=utils.Hlt1Lines_Wmumuqq(), Hlt2_decisions_=utils.Hlt2Lines_Wmumuqq())
    mc_truth_vars = FunctorCollection({})
    mctruth = None
    if options.simulation == True:
        #all_vars += utils.make_MCTruth_vars(data, algname)
        mc_truth_vars_temp, mctruth = utils.make_MCTruth_vars(data, algname)
        mc_truth_vars += mc_truth_vars_temp

    track_vars = utils.make_track_vars(pvs)


    variables = {
        "Mu": all_vars+track_vars+mc_truth_vars+utils.make_isolation_vars(data, "Mu", "W"),
    }

    evt_vars = utils.make_evt_vars(pvs,  Hlt1_decisions_=utils.Hlt1Lines_Wmumuqq(), Hlt2_decisions_=utils.Hlt2Lines_Wmumuqq(), is_simulation=options.simulation, Spruce_decisions_ = utils.SpruceLines_Wmumuqq())

    funtuple = Funtuple(
        name=algname,
        tuple_name="DecayTree",
        fields=fields,
        variables=variables,
        event_variables=evt_vars,
        inputs=data,
    )

    mct_fields = {
        "Mu" :   "[mu-]CC",
    }


    mc_particles = get_mc_particles('/Event/MC/Particles')

    variables_mc = utils.make_mc_variables(mc_particles)

    mct_variables = { 'ALL': variables_mc}

    funtuple_mc = FuntupleMC(
        name=algname+"MC",
        tuple_name="DecayTree",
        fields=mct_fields,
        variables=mct_variables,
        inputs=mc_particles,
    )

    if options.simulation == True:
        algs = {algname+"MC": [require_pvs(pvs), funtuple, funtuple_mc]}
        return algs
    else:
        algs = {algname: [require_pvs(pvs), funtuple]}
        return algs

def main(options: Options):

    name = "myTupleWmunu"

    algs = make_algs_Wmunu(options, name)

    return make_config(options, algs)
