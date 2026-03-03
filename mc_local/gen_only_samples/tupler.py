from FunTuple import FunctorCollection, FunTuple_MCParticles as FuntupleMC
import FunTuple.functorcollections as FC
from PyConf.Algorithms import PrintMCTree
import Functors as F
from Functors.math import in_range
from DaVinci import Options, make_config
from PyConf.reading import get_mc_particles, get_mc_header, get_pvs




def main(options: Options):
    # FunTuple: define fields (branches)
    fields = {
        "Bcp": "[B_c+ -> (H_10 -> pi+ mu-) mu+]CC",
        "HNL": "[B_c+ -> ^(H_10 -> pi+ mu-) mu+]CC",
        "PipHNL": "[B_c+ -> (H_10 -> ^pi+ mu-) mu+]CC",
        "MumHNL": "[B_c+ -> (H_10 -> pi+ ^mu-) mu+]CC",
        "Mup": "[B_c+ -> (H_10 -> pi+ mu-) ^mu+]CC",
    }

    pvs = get_pvs()

    # FunTuple: define common variables
    variables_all = FunctorCollection({
        "PT": F.PT, 
        "ETA": F.ETA, 
        "PHI": F.PHI, 
        "P": F.P, 
        "E" : F.ENERGY,
        "M": F.MASS,
        })

    vertex_vars = FunctorCollection({
        "END_VX": F.END_VX,
        "END_VY": F.END_VY,
        "END_VZ": F.END_VZ,
        "END_VRHO": F.END_VRHO,
    })

    # FunTuple: associate functor collections to field (branch) name
    variables = {
        "ALL": variables_all,
        "Bcp": vertex_vars,
        "HNL": vertex_vars,
    }

    # FunTuple: define input data
    all_mc_particles = get_mc_particles("/Event/MC/Particles")

    # FunTuple: define event-level variables using functor collections
    mc_header = get_mc_header()
    evt_vars = FC.MCPrimaries(mc_header=mc_header)

    printMC = PrintMCTree(
        MCParticles=all_mc_particles, ParticleNames=["B_c+", "B_c-"], OutputLevel=4
    )

    # now we filter that set of B particles
    cut = F.require_all(in_range(1.8, F.ETA, 5.0))

    tuple_1 = FuntupleMC(
        name="myTuple",
        tuple_name="DecayTree",
        fields=fields,
        variables=variables,
        event_variables=evt_vars,
        inputs=all_mc_particles,
    )

    return make_config(options, {"Tuple": [printMC, tuple_1]})
