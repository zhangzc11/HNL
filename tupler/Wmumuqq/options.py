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


def add_nJ_text(text, nJ):
    if nJ <= 0:
        return text
    replacement = " ".join(["CELLjet"] * nJ)
    return text.replace("CELLjet", replacement)

def add_nJ(fields, nJ):
    result = {}
    for key, value in fields.items():
        result[key] = add_nJ_text(value, nJ)
    return result

def mark_ith_CELLjet(text, i, target="CELLjet", replacement="^CELLjet"):
    if i < 1:
        return text
    parts = text.split(target)
    if len(parts) <= i:
        return text
    result = target.join(parts[:i]) + replacement + target.join(parts[i:i+1])
    if i < len(parts) - 1:
        result += target + target.join(parts[i+1:])
    return result

def convert_decay(fields, isSS, isMCT, nJ):
    fields = add_nJ(fields, nJ)
    base_text = fields["W"]
    for iJ in range(nJ):
        fields["Jet"+str(iJ+1)] = mark_ith_CELLjet(base_text, iJ+1, "CELLjet", "^CELLjet")

    #print(fields)
    result = {}
    for key, value in fields.items():
        new_value = value
        if isSS:
            new_value = new_value.replace("mu+", "mu-")
        if isMCT:
            new_value = new_value.replace("->", "=>").replace("CELLjet", "Hadron").replace("W-", "W_R-")
        result[key] = new_value
    return result

def decay_descriptor(isSS = False, isMCT = False, nJ=1):
    fields = {
        "W" :   "[W- -> mu- (nu_Rmu -> mu+ CELLjet)]CC",
        "Mu" :   "[W- -> ^mu- (nu_Rmu -> mu+ CELLjet)]CC",
        "NuR" :   "[W- -> mu- ^(nu_Rmu -> mu+ CELLjet)]CC",
        "MuNuR" :   "[W- -> mu- (nu_Rmu -> ^mu+ CELLjet)]CC",
    }
    return convert_decay(fields, isSS, isMCT, nJ)

def diyline(isSS=False, nJ=1):
    name_suffix = "_OS_"
    if isSS:
        name_suffix = "_SS_"
    name_suffix = name_suffix + str(nJ)+"J"

    muons_W = ParticleFilter(
        standard_particles.make_long_muons(),
        F.FILTER(F.require_all(F.PT > 15 * GeV, F.ISMUON == True, F.INECAL == True, F.INHCAL == True))
        )

    muons_long = standard_particles.make_long_muons()
    muons_down = standard_particles.make_down_muons()
    #muons_Ttrack = standard_particles.make_ttrack_muons(make_protoparticles=utils.make_ttrack_protoparticles_diy)
    
    #muons_all = ParticleContainersMerger([muons_long, muons_down, muons_Ttrack], name="muons_N_all")
    muons_all = ParticleContainersMerger([muons_long, muons_down], name="muons_N_all")

    muons_N = ParticleFilter(
        muons_all,
        F.FILTER(F.require_all(F.PT > 2 * GeV, F.ISMUON == True))
        )

    jets = utils.make_jets(pt_min = utils.JET_PT_MIN)

    combination_code_HNL = F.require_all(
        F.MASS > 0.1 * GeV,
        F.PT > 10. * GeV,
    )

    vertex_code_HNL = F.require_all(
        F.MASS > 0.1 * GeV,
        F.CHARGE > -99,
    )

    decay_HNL = '[nu_Rmu -> mu+ CELLjet]cc'
    if isSS:
        decay_HNL = '[nu_Rmu -> mu- CELLjet]cc'
    decay_HNL = add_nJ_text(decay_HNL, nJ)
    input_HNL = [muons_N]
    for iJ in range(nJ):
        input_HNL.append(jets)

    HNL = ParticleCombiner(
        name="myHNL"+name_suffix,
        Inputs=input_HNL,
        DecayDescriptor=decay_HNL,
        CombinationCut=combination_code_HNL,
        CompositeCut=vertex_code_HNL)
 
    combination_code_W = F.require_all(
        fmath.in_range( 20.0 * GeV, F.MASS, 200.0 * GeV),
        utils.make_comb_mass(1) > 10.0 * GeV,
    )

    vertex_code_W = F.require_all(
        F.MASS > 0.1 * GeV,
        F.CHARGE > -99,
    )
       
    WR = ParticleCombiner(
        name="myWR"+name_suffix,
        Inputs=[muons_W, HNL],
        DecayDescriptor='[W- -> mu- nu_Rmu]cc',
        CombinationCut=combination_code_W,
        CompositeCut=vertex_code_W)  

    return WR, muons_W, muons_N
  

def make_algs_Wmumuqq(options, algname, isSS_=False, nJ_=1, doReco=True, mcInclusive=False):
    data, muons_W, muons_N = diyline(isSS_, nJ_)

    pvs = get_pvs()

    fields = decay_descriptor(isSS = isSS_, isMCT = False, nJ=nJ_)

    all_vars = utils.make_base_vars(pvs)
    all_vars += utils.make_HltTisTos_vars(data, Hlt1_decisions_=utils.Hlt1Lines_Wmumuqq(), Hlt2_decisions_=utils.Hlt2Lines_Wmumuqq())
    mc_truth_vars = FunctorCollection({})
    mctruth = None
    if options.simulation == True and mcInclusive == False:
        #all_vars += utils.make_MCTruth_vars(data, algname)
        mc_truth_vars_temp, mctruth = utils.make_MCTruth_vars(data, algname)
        mc_truth_vars += mc_truth_vars_temp

    vertex_vars = utils.make_vertex_vars(pvs)
    vertex_vars_3daughters = utils.make_vertex_vars_3daughters()
    vertex_vars_4daughters = utils.make_vertex_vars_4daughters()
    track_vars = utils.make_track_vars(pvs)

    jet_vars = utils.make_jet_vars(data, mctruth, is_simulation=options.simulation)
    wcomb_vars = utils.make_WR_comb_vars()

    variables = {
        #"ALL": all_vars,
        "W": all_vars+vertex_vars+mc_truth_vars+wcomb_vars,
        "NuR": all_vars+vertex_vars+mc_truth_vars,
        "Mu": all_vars+track_vars+mc_truth_vars+utils.make_isolation_vars(muons_W, "Mu", "W"),
        "MuNuR": all_vars+track_vars+mc_truth_vars+utils.make_isolation_vars(muons_N, "Mu", "W"),
    }
    for iJ in range(nJ_):
        variables["Jet"+str(iJ+1)] = jet_vars

    evt_vars = utils.make_evt_vars(pvs,  Hlt1_decisions_=utils.Hlt1Lines_Wmumuqq(), Hlt2_decisions_=utils.Hlt2Lines_Wmumuqq(), is_simulation=options.simulation, Spruce_decisions_ = utils.SpruceLines_Wmumuqq())

    funtuple = Funtuple(
        name=algname,
        tuple_name="DecayTree",
        fields=fields,
        variables=variables,
        event_variables=evt_vars,
        inputs=data,
    )

    mct_fields = decay_descriptor(isSS = isSS_, isMCT = True, nJ=nJ_)
    if mcInclusive:
        mu_ = "mu+"
        if isSS_:
            mu_ = "mu-"
        mct_fields = {
            #"W": "[W_R- => mu- nu_Rmu]CC", # 999/1000 
            #"Mu": "[W_R- => ^mu- nu_Rmu]CC", # 999/1000 
            #"NuR": "[W_R- => mu- ^nu_Rmu]CC", # 999/1000 
            #"W": "[W_R- => mu- (nu_Rmu -> X ... ) ]CC", # 999/1000
            #"W": "[W_R- => mu- (nu_Rmu -> Lepton ... ) ]CC", # 640/1000
            #"W": "[W_R- => mu- (nu_Rmu -> mu- ... ) ]CC", # 318/1000
            #"W": "[W_R- => mu- (nu_Rmu -> mu+ ... ) ]CC", # 222/1000
            #"W": "[W_R- => mu- (nu_Rmu -> X X ... ) ]CC", # 640 /1000
            #"W": "[W_R- => mu- (nu_Rmu => X) ]CC", # 359/1000, X is nu_Rmu, not decayed?
            #"NuR": "[W_R- => mu- ^(nu_Rmu => X) ]CC", # 359/1000, X is nu_Rmu, not decayed?
            #"X": "[W_R- => mu- (nu_Rmu => ^X) ]CC", # 359/1000, X is nu_Rmu, not decayed?
            #"W": "[W_R- => mu- (nu_Rmu => X X) ]CC", # 34/1000
            #"W": "[W_R- => mu- (nu_Rmu => X X X) ]CC", # 215/1000
            #"W": "[W_R- => mu- (nu_Rmu => X X X X) ]CC", # 276/1000
            #"W": "[W_R- => mu- (nu_Rmu => X X X X X) ]CC", # 104/1000
            #"W": "[W_R- => mu- (nu_Rmu => X X X X X X) ]CC", # 44/1000
            #"W": "[W_R- => mu- (nu_Rmu => Lepton Hadron ) ]CC", # 34/1000
            "W": "[W_R- => mu- (nu_Rmu -> "+mu_+" Hadron ... ) ]CC", # 637/1000
            "Mu": "[W_R- => ^mu- (nu_Rmu -> "+mu_+" Hadron ... ) ]CC", # 637/1000
            "NuR": "[W_R- => mu- ^(nu_Rmu -> "+mu_+" Hadron ... ) ]CC", # 637/1000
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
        algs = {algname+"MC": [funtuple_mc]}
        if doReco:
            algs.update({algname: [require_pvs(pvs), funtuple]})
        return algs
    else:
        algs = {algname: [require_pvs(pvs), funtuple]}
        return algs

def main(options: Options):

    name = "myTuple"

    algs = make_algs_Wmumuqq(options, name+"SS1J", isSS_=True, nJ_=1)
    algs.update(make_algs_Wmumuqq(options, name+"SS2J", isSS_=True, nJ_=2))
    algs.update(make_algs_Wmumuqq(options, name+"OS1J", isSS_=False, nJ_=1))
    algs.update(make_algs_Wmumuqq(options, name+"OS2J", isSS_=False, nJ_=2))
    if options.simulation == True:
        algs.update(make_algs_Wmumuqq(options, name+"SS3J", isSS_=True, nJ_=3, doReco=False))
        algs.update(make_algs_Wmumuqq(options, name+"SS4J", isSS_=True, nJ_=4, doReco=False))
        algs.update(make_algs_Wmumuqq(options, name+"SS5J", isSS_=True, nJ_=5, doReco=False))
        algs.update(make_algs_Wmumuqq(options, name+"OS3J", isSS_=False, nJ_=3, doReco=False))
        algs.update(make_algs_Wmumuqq(options, name+"OS4J", isSS_=False, nJ_=4, doReco=False))
        algs.update(make_algs_Wmumuqq(options, name+"OS5J", isSS_=False, nJ_=5, doReco=False))
        algs.update(make_algs_Wmumuqq(options, name+"SSInclusive", isSS_=True, nJ_=1, doReco=False, mcInclusive=True))
        algs.update(make_algs_Wmumuqq(options, name+"OSInclusive", isSS_=False, nJ_=1, doReco=False, mcInclusive=True))

    return make_config(options, algs)
