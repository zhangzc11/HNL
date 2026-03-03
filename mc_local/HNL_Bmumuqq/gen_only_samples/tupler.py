from FunTuple import FunctorCollection, FunTuple_MCParticles as FuntupleMC
import FunTuple.functorcollections as FC
from PyConf.Algorithms import PrintMCTree
import Functors as F
from Functors.math import in_range
import Functors.math as fmath
from DaVinci import Options, make_config
from PyConf.reading import get_mc_particles, get_mc_header, get_pvs




def main(options: Options):
    # FunTuple: define fields (branches)
    fields = {
        "B0": "[B0 => mu+ X (H_10 => Lepton X)]CC",
        "Pi": "[B0 => mu+ ^X (H_10 => Lepton X)]CC",
        "Mu": "[B0 => ^mu+ X (H_10 => Lepton X)]CC",
        "NuR": "[B0 => mu+ X ^(H_10 => Lepton X)]CC",
        "MuNuR": "[B0 => mu+ X (H_10 => ^Lepton X)]CC",
        "Jet1": "[B0 => mu+ X (H_10 => Lepton ^X)]CC",
    }

    fields1 = {"B0": "B0 => mu+ X (H_10 => mu+ X)"}
    #fields1 = {"B0": "B0 => mu+ X (H_10 -> mu+ X ... )"}
    fields2 = {"B0": "B0 => mu+ X (H_10 => mu- X)"}
    fields3 = {"B0": "B~0 => mu- X (H_10 => mu+ X)"}
    fields4 = {"B0": "B~0 => mu- X (H_10 => mu- X)"}

    pvs = get_pvs()

    # FunTuple: define common variables
    variables_all = FunctorCollection({
        "PT": F.PT, 
        "ETA": F.ETA, 
        "PHI": F.PHI, 
        "P": F.P, 
        "E" : F.ENERGY,
        "M": F.MASS,
        "ID": F.PARTICLE_ID,
        "Q": F.CHARGE,
        "THETA": F.THETA,
        })

    vertex_vars = FunctorCollection({
        "END_VX": F.END_VX,
        "END_VY": F.END_VY,
        "END_VZ": F.END_VZ,
        "END_VR": fmath.sqrt(F.END_VZ*F.END_VZ + F.END_VY*F.END_VY + F.END_VX*F.END_VX),
        "END_VRHO": F.END_VRHO,
        "ORIGIN_VX": F.ORIGIN_VX,
        "ORIGIN_VY": F.ORIGIN_VY,
        "ORIGIN_VZ": F.ORIGIN_VZ,
    })

    # FunTuple: associate functor collections to field (branch) name
    variables = {
        "ALL": variables_all+vertex_vars,
    }

    # FunTuple: define input data
    all_mc_particles = get_mc_particles("/Event/MC/Particles")

    # FunTuple: define event-level variables using functor collections
    mc_header = get_mc_header()
    evt_vars = FC.MCPrimaries(mc_header=mc_header)

    printMC = PrintMCTree(
        MCParticles=all_mc_particles, ParticleNames=["B0", "B0"], OutputLevel=4
    )

    # now we filter that set of B particles
    cut = F.require_all(in_range(1.8, F.ETA, 5.0))

    tuple = FuntupleMC(
        name="myTuple",
        tuple_name="DecayTree",
        fields=fields,
        variables=variables,
        event_variables=evt_vars,
        inputs=all_mc_particles,
    )

    tuple_1 = FuntupleMC(
        name="myTuple1",
        tuple_name="DecayTree",
        fields=fields1,
        variables=variables,
        event_variables=evt_vars,
        inputs=all_mc_particles,
    )

    tuple_2 = FuntupleMC(
        name="myTuple2",
        tuple_name="DecayTree",
        fields=fields2,
        variables=variables,
        event_variables=evt_vars,
        inputs=all_mc_particles,
    )

    tuple_3 = FuntupleMC(
        name="myTuple3",
        tuple_name="DecayTree",
        fields=fields3,
        variables=variables,
        event_variables=evt_vars,
        inputs=all_mc_particles,
    )

    tuple_4 = FuntupleMC(
        name="myTuple4",
        tuple_name="DecayTree",
        fields=fields4,
        variables=variables,
        event_variables=evt_vars,
        inputs=all_mc_particles,
    )

    return make_config(options, {"Tuple": [tuple], "Tuple1": [tuple_1], "Tuple2": [tuple_2], "Tuple3": [tuple_3], "Tuple4": [tuple_4]})
