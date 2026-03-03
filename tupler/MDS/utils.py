from FunTuple import FunctorCollection
import FunTuple.functorcollections as FC
import Functors as F
from DecayTreeFitter import DecayTreeFitter
from DaVinciMCTools import MCTruthAndBkgCat
from PyConf.reading import get_rec_summary, get_odin
import Functors.math as fmath
from RecoConf.muonid import make_muon_hits
from Hlt2Conf.lines.rd.builders.rd_isolation import find_in_decay
from PyConf.Algorithms import LHCb__JetInfoRelationTable as JetInfoRelationTable
from PyConf.reading import get_mc_track_info
from DaVinciMCTools import MCReconstructed as MCRected
from DaVinciMCTools import MCReconstructible as MCRectible

from PyConf.Algorithms import WeightedRelTableAlg
from Hlt2Conf.standard_jets import make_particleflow
from RecoConf.algorithms_thor import ParticleFilter
from GaudiKernel.SystemOfUnits import GeV
from Hlt2Conf.standard_jets import build_jets
from PyConf.Algorithms import ParticleFlowMakerMC
from RecoConf.data_from_file import mc_unpackers
from PyConf.reading import get_mc_particles

JET_PT_MIN = 2.0 * GeV

def make_base_vars(pvs):
    base_vars = FunctorCollection({
        "Key": F.OBJECT_KEY,
        "M": F.MASS,
        "P": F.P,
        "log_P": fmath.log(F.P),
        "E": F.ENERGY,
        "ETA": F.ETA,
        "PHI": F.PHI,
        "THETA": F.THETA,
        "PT": F.PT,
        "log_PT": fmath.log(F.PT),
        "PX": F.PX,
        "PY": F.PY,
        "PZ": F.PZ,
        "ID": F.PARTICLE_ID,
        "Q": F.CHARGE,
        "CHI2": F.CHI2,                 # chi2
        "CHI2DOF": F.CHI2DOF,           # chi2 of degrees of freedom
        "OWNPVIP": F.OWNPVIP,           # Impact parameter wrt own PV
        "OWNPVIPCHI2": F.OWNPVIPCHI2,   # Impact parameter chi2 wrt own PV
        "OWNPV_X": F.OWNPVX,            # x-coordinate of best PV
        "OWNPV_Y": F.OWNPVY,            # y-coordinate of best PV
        "OWNPV_Z": F.OWNPVZ,
        "OWNPVNDOF": F.VALUE_OR(-1) @ F.NDOF @ F.OWNPV,
        "MINIPCHI2": F.MINIPCHI2(pvs),
        "MINIP": F.MINIP(pvs),
        "ALLPV_IP[nPVs]": F.ALLPV_IP(pvs),
        "ALLPV_IPCHI2[nPVs]": F.MAP(F.IP).bind(F.MAP(F.TOLINALG @ F.POSITION) @ F.TES(pvs), F.FORWARDARGS),
        })
    base_vars += FC.ParticleID()
    base_vars += FC.NeutralCaloInfo(extra_info=True)
    base_vars += FC.ChargedCaloInfo(extra_info=True)

    return base_vars

def make_vertex_vars(pvs):
    return FunctorCollection({
        "VTXCHI2NDOF": F.CHI2DOF,         # Vertex fit chi2/ndf
        "END_VX": F.END_VX,               # x-coordinate of decay vertex
        "END_VY": F.END_VY,               # y-coordinate of decay vertex
        "END_VZ": F.END_VZ,               # z-coordinate of decay vertex
        "BPV_LTIME": F.BPVLTIME(pvs),      # Lifetime with respect to the BPV
        "BPV_IPCHI2": F.BPVIPCHI2(pvs),    # Impact parameter Chi2 with respect to the BPV
        # OWNPV values
        "OWNPV_DIRA": F.OWNPVDIRA,        # Direction angle cosine wrt own PV
        "OWNPV_FD": F.OWNPVFD,            # Flight distance wrt own PV
        "OWNPV_FDCHI2": F.OWNPVFDCHI2,    # Flight distance chi2 wrt own PV
        #"OWNPV_IPCHI2": F.OWNPVIPCHI2,    # Impact parameter Chi2 with respect to the own PV
        "OWNPV_FDIR": F.OWNPVFDIR,         # unity vector of OWNPVFDVEC
        "OWNPV_FDVEC": F.OWNPVFDVEC,       # Three-vector distance between the endvertex position and the position of the PV associated to a particle
        "OWNPV_VDRHO": F.OWNPVVDRHO,      # Radial flight distance wrt own PV
        "OWNPV_VDX": F.OWNPVVDX,
        "OWNPV_VDY": F.OWNPVVDY,
        "OWNPV_VDZ": F.OWNPVVDZ,          # z-direction flight distance
        "OWNPV_LTIME": F.OWNPVLTIME,      # Proper lifetime
        "OWNPV_DLS": F.OWNPVDLS,          # Decay length significance
        "OWNPV_ETA": F.OWNPVETA,
        "ALLPV_FD[nPVs]": F.ALLPV_FD(pvs),
        "ALV": F.ALV(Child1=1, Child2=2),
        # DOCA
        "DOCA12": F.DOCA(1, 2),           # DOCA between first and second daughter
        "SDOCA12": F.SDOCA(1, 2),         # SDOCA between first and second daughter
        "DOCA12CHI2": F.DOCACHI2(1, 2),   # DOCA chi2 between first and second daughter
        # Daughter Max, Min and Sums
        "MAX_PT": F.MAX(F.PT),            # Maximum PT of daughters
        "log_MAX_PT": fmath.log(F.MAX(F.PT)),            # Maximum PT of daughters
        "MIN_PT": F.MIN(F.PT),            # Minimum PT of daughters
        "log_MIN_PT": fmath.log(F.MIN(F.PT)),            # Minimum PT of daughters
        "SUM_PT": F.SUM(F.PT),            # Sum of daughters' PT
        "log_SUM_PT": fmath.log(F.SUM(F.PT)),            # Sum of daughters' PT
        "MAX_P": F.MAX(F.P),              # Maximum momentum of daughters
        "log_MAX_P": fmath.log(F.MAX(F.P)),              # Maximum momentum of daughters
        "MIN_P": F.MIN(F.P),              # Minimum momentum of daughters
        "log_MIN_P": fmath.log(F.MIN(F.P)),              # Minimum momentum of daughters
        "SUM_P": F.SUM(F.P),              # Sum of daughters' momentum
        "log_SUM_P": fmath.log(F.SUM(F.P)),              # Sum of daughters' momentum
        "MAX_OWNPVIPCHI2": F.MAX(F.OWNPVIPCHI2),  # Max IP chi2 of daughters
        "MIN_OWNPVIPCHI2": F.MIN(F.OWNPVIPCHI2),  # Min IP chi2 of daughters
        "SUM_OWNPVIPCHI2": F.SUM(F.OWNPVIPCHI2),  # Sum of daughters' IP chi2
        "MAX_DOCACHI2": F.MAXDOCACHI2,      # Maximum DOCA chi2 between any daughters
        "MAX_DOCA": F.MAXDOCA,              # Maximum DOCA between any daughters
        "MAX_SDOCACHI2": F.MAXSDOCACHI2,    # Maximum signed DOCA chi2
        "MAX_SDOCA": F.MAXSDOCA,            # Maximum signed DOCA
        "MAX_MINIPCHI2": F.MAX(F.MINIPCHI2(pvs)),
        "SUM_MINIPCHI2": F.SUM(F.MINIPCHI2(pvs)),
        "OBJECT_KEY": F.OBJECT_KEY,
        "KEY": F.VALUE_OR(-1) @ F.OBJECT_KEY @ F.TRACK,
    })

def make_vertex_vars_3daughters():
    return FunctorCollection({
        "DOCA13": F.DOCA(1, 3),
        "SDOCA13": F.SDOCA(1, 3),
        "DOCA13CHI2": F.DOCACHI2(1, 3),
        "DOCA23": F.DOCA(2, 3),
        "SDOCA23": F.SDOCA(2, 3),
        "DOCA23CHI2": F.DOCACHI2(2, 3),
    })

def make_vertex_vars_4daughters():
    return FunctorCollection({
        "DOCA13": F.DOCA(1, 3),
        "SDOCA13": F.SDOCA(1, 3),                                                                                                   "DOCA13CHI2": F.DOCACHI2(1, 3),
        "DOCA14": F.DOCA(1, 4),
        "SDOCA14": F.SDOCA(1, 4),
        "DOCA14CHI2": F.DOCACHI2(1, 4),
        "DOCA23": F.DOCA(2, 3),
        "SDOCA23": F.SDOCA(2, 3),
        "DOCA23CHI2": F.DOCACHI2(2, 3),
        "DOCA24": F.DOCA(2, 4),
        "SDOCA24": F.SDOCA(2, 4),
        "DOCA24CHI2": F.DOCACHI2(2, 4),
        "DOCA34": F.DOCA(3, 4),
        "SDOCA34": F.SDOCA(3, 4),
        "DOCA34CHI2": F.DOCACHI2(3, 4),
    })


def make_track_vars(pvs):
    return  FunctorCollection({
        "TRACKPT": F.TRACK_PT,
        # Standard PID
        "PIDp": F.PID_P,              # Proton PID likelihood
        "PIDK": F.PID_K,              # Kaon PID likelihood
        "PIDPi": F.PID_PI,            # Pion PID likelihood
        "PIDe": F.PID_E,              # Electron PID likelihood
        "PIDmu": F.PID_MU,            # Muon PID likelihood
        "GHOSTPROB": F.GHOSTPROB,
        "ISMUON": F.ISMUON,              # Boolean: is it identified as a muon 0 or 1?
        "INMUON": F.INMUON,
        "TRCHI2":    F.CHI2 @ F.TRACK,
        "TRCHI2DOF": F.CHI2DOF @ F.TRACK,
        "QOVERP": F.QOVERP @ F.TRACK,
        "TRNDOF": F.VALUE_OR(-1) @ F.NDOF @ F.TRACK,                # NDOF in track fit, if this is higher then more hits used in fit
        "NHITS": F.VALUE_OR(-1) @ F.NHITS @ F.TRACK,                # Total number of hits in all detectors
        "NFTHITS": F.VALUE_OR(-1) @ F.NFTHITS @ F.TRACK,            # Total number of hits in FT (SciFi)
        "NVPHITS": F.VALUE_OR(-1) @ F.NVPHITS @ F.TRACK,            # Total number of hits in VELO phi sensors
        "NVPHITSC":  F.VALUE_OR(-1) @ F.NVPHITSC @ F.TRACK,         # Total number of hits in VELO C-side
        "NVPHITSA":  F.VALUE_OR(-1) @ F.NVPHITSA @ F.TRACK,         # Total number of hits in VELO A-side
        "NVPHOLES":  F.VALUE_OR(-1) @ F.NVPHOLES @ F.TRACK,         # Total number of holes in VELO of track as defined by HitPattern
        "NFTHOLES":  F.VALUE_OR(-1) @ F.NFTHOLES @ F.TRACK,         # Total number of holes in SciFi of track as defined by HitPattern
        "NVPLAYERS":  F.VALUE_OR(-1) @ F.NVPLAYERS @ F.TRACK,       # Total number of VELO layers of track as defined by HitPattern
        "NUTHITS": F.VALUE_OR(-1) @ F.NUTHITS @ F.TRACK,            # Total number of hits in UT
        "NFTHITS": F.VALUE_OR(-1) @ F.NFTHITS @ F.TRACK,            # Total number of hits in Fibre Tracker (SciFi)
        "TRACKHISTORY": F.VALUE_OR(-1) @ F.TRACKHISTORY @ F.TRACK,  # Track reconstruction history
        "TRACKHASVELO": F.VALUE_OR(-1) @ F.TRACKHASVELO @ F.TRACK,
        "TRACKHASUT": F.VALUE_OR(-1) @ F.TRACKHASUT @ F.TRACK,
        "TRACKHAST": F.VALUE_OR(-1) @ F.TRACKHAST @ F.TRACK,
        "PPHASRICH": F.VALUE_OR(-1) @ F.PPHASRICH @ F.PROTOPARTICLE,
        "ISLONG": F.VALUE_OR(-1) @ F.TRACKISLONG @ F.TRACK,
        "ISDOWNSTREAM": F.VALUE_OR(-1) @ F.TRACKISDOWNSTREAM @ F.TRACK,
        "ISTTRACK": F.VALUE_OR(-1) @ F.TRACKISTTRACK @ F.TRACK,
        "ISUPSTREAM": F.VALUE_OR(-1) @ F.TRACKISUPSTREAM @ F.TRACK,
        "ISVELO": F.VALUE_OR(-1) @ F.TRACKISVELO @ F.TRACK,
        "ISCLONE": F.VALUE_OR(-1) @ F.TRACKISCLONE @ F.TRACK,
        "TX": F.TX,
        "TY": F.TY,
        "POSITION_STATEAT_FirstMeasurement_X"     : F.POSITION_X @ F.STATE_AT("FirstMeasurement")@ F.TRACK,
        "POSITION_STATEAT_FirstMeasurement_Y"     : F.POSITION_Y @ F.STATE_AT("FirstMeasurement")@ F.TRACK,
        "POSITION_STATEAT_FirstMeasurement_Z"     : F.POSITION_Z @ F.STATE_AT("FirstMeasurement")@ F.TRACK,
        "POSITION_STATEAT_LastMeasurement_X"      : F.POSITION_X @ F.STATE_AT("LastMeasurement") @ F.TRACK,
        "POSITION_STATEAT_LastMeasurement_Y"      : F.POSITION_Y @ F.STATE_AT("LastMeasurement") @ F.TRACK,
        "POSITION_STATEAT_LastMeasurement_Z"      : F.POSITION_Z @ F.STATE_AT("LastMeasurement") @ F.TRACK,
        })


def make_DTF_vars(data_, mass_to_constrain_, name_="DTF", name2_=""):
    dtf = DecayTreeFitter(name=name_+"_"+name2_, input_particles=data_, mass_constraints=mass_to_constrain_) 
    dtf_ownpv = DecayTreeFitter(name=name_+"_"+name2_+"_OWNPV", input_particles=data_, mass_constraints=mass_to_constrain_, constrain_to_ownpv=True) 
    dtf_vars_base = FunctorCollection({
        "CHI2DOF" : F.CHI2DOF,
        "OWNPV_DIRA": F.OWNPVDIRA,        # Direction angle cosine wrt own PV
        "OWNPV_FD": F.OWNPVFD,            # Flight distance wrt own PV
        "OWNPV_FDCHI2": F.OWNPVFDCHI2,    # Flight distance chi2 wrt own PV
        "OWNPVIP": F.OWNPVIP,           # Impact parameter wrt own PV
        "OWNPVIPCHI2": F.OWNPVIPCHI2,   # Impact parameter chi2 wrt own PV
        "ETA" : F.ETA,
        "PHI" : F.PHI,
        "M" : F.MASS,
        "P" : F.P,
        "PT" : F.PT,
        "PX" : F.PX,
        "PY" : F.PY,
        "PZ" : F.PZ,
        "E" : F.ENERGY,
        })
    dtf_vars = FunctorCollection(
       {name_+"_" + k: dtf(v) for k, v in dtf_vars_base.get_thor_functors().items()}
    )

    dtf_vars += FunctorCollection(
       {name_+"_OWNPV_" + k: dtf_ownpv(v) for k, v in dtf_vars_base.get_thor_functors().items()},
    )

    return dtf_vars



def make_particleflowMC(
        requestedParticlesPIDs=[11, -11, 13, -13, 211, -211, 321, -321, 2112, -2112, 2212, -2212, 22, 310, 311, 130, 111, 3122, -3122], 
        banPIDs=[12, -12, 14, -14, 16, -16, 23]): # ν_e, ν̄_e, ν_μ, ν̄_μ, ν_τ, ν̄_τ, Z0
    return ParticleFlowMakerMC(
        Inputs=[mc_unpackers()["MCParticles"]],
        requestedParticlesPIDs=requestedParticlesPIDs,
        banPIDs=banPIDs).Output

def make_jets(name="SimpleJets_{hash}", pt_min=10 * GeV, JetsByVtx=False, useMC=False):
    """
    Create jets from the particle flow, optionally using MC information.
    Note: JetsByVtx  = True enforces jets created per vertex in the event, switching to False enforces jet production with all particles in an event.
    """
    if useMC:
        pflow = make_particleflowMC()
        jets = build_jets(pflow=pflow, MCJets=useMC, JetsByVtx=JetsByVtx, name="MCJetBuilder" + name)
    else:
        pflow = make_particleflow()
        jets = build_jets(pflow=pflow, JetsByVtx=JetsByVtx, MCJets=useMC, name="JetBuilder" + name)

    code = F.require_all(
        F.IS_ABS_ID("CELLjet"), F.PT > pt_min, F.NINGENERATION(F.CHARGE != 0, 1) > 1
    )
    return ParticleFilter(jets, F.FILTER(code), name=name)

def matchjets(refjet, inputjet):
    """
    Match MC jets to reco jets using a DR cut.
    """

    alg = WeightedRelTableAlg(
        name="MatchMCjets_{hash}",
        ReferenceParticles=refjet,
        InputCandidates=inputjet,
        Cut=F.require_all(F.DR2 < (0.5**2)),
    )
    ETA_EXTRA = (F.ETA_COORDINATE @ F.SLOPES) @ F.TO @ F.FORWARDARG0
    PHI_EXTRA = (F.PHI_COORDINATE @ F.SLOPES) @ F.TO @ F.FORWARDARG0
    ETA_REF = (F.ETA_COORDINATE @ F.SLOPES)  @ F.FORWARDARG1
    PHI_REF = (F.PHI_COORDINATE @ F.SLOPES) @ F.FORWARDARG1
    DETA = (ETA_REF - ETA_EXTRA)
    DPHI = F.ADJUST_ANGLE  @ (PHI_REF - PHI_EXTRA)

    DR2 = DETA * DETA + DPHI * DPHI

    functors = FunctorCollection(
        {
        "MATCH": F.VALUE_OR(0) @ F.MAP_INPUT_SIZE(alg.OutputRelations),
        #"MCJets_PT":  F.VALUE_OR(-99) @ F.MAP_INPUT_ARRAY(Functor=F.PT, Relations=alg.OutputRelations),
        "MCClosestJet_PT":  F.VALUE_OR(-99) @ F.PT @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_ETA":  F.VALUE_OR(-99) @ F.ETA @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_PHI":  F.VALUE_OR(-99) @ F.PHI @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_ID":  F.VALUE_OR(-99) @ F.PARTICLE_ID @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        }
    )

    return functors

def JetConstituents(data__, MCTRUTH):
    """
    Get the constituents of the jets.
    """
    #MCTRUTH = MCTruthAndBkgCat(input_particles=data__, name="MCTruthAndBkgCat_functor_jets_{hash}")
    #MCMOTHER_ID = lambda n: F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(n, F.PARTICLE_ID))

    return FunctorCollection({
        "nConstituents": F.NINGENERATION(F.ALL,1),
        "nChargedConstituents": F.NINGENERATION(F.CHARGE != 0, 1),
        # Store an array of jet constituents' Px, Py, Pz
        "ConstituentPx": F.MAP(F.PX) @ F.GET_GENERATION(1),
        "ConstituentPy": F.MAP(F.PY) @ F.GET_GENERATION(1),
        "ConstituentPz": F.MAP(F.PZ) @ F.GET_GENERATION(1),
        "ConstituentPt": F.MAP(F.PT) @ F.GET_GENERATION(1),
        "ConstituentEta": F.MAP(F.ETA) @ F.GET_GENERATION(1),
        "ConstituentPhi": F.MAP(F.PHI) @ F.GET_GENERATION(1),
        "ConstituentMass": F.MAP(F.MASS) @ F.GET_GENERATION(1),
        "ConstituentPID": F.MAP(F.PARTICLE_ID) @ F.GET_GENERATION(1),
        "ConstituentTRUEID": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.PARTICLE_ID) ) @ F.GET_GENERATION(1),
        "ConstituentMC_MOTHER_ID": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(1, F.PARTICLE_ID)) ) @ F.GET_GENERATION(1),
        "ConstituentMC_MOTHER_Key": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(1, F.OBJECT_KEY)) ) @ F.GET_GENERATION(1),
        "ConstituentMC_GD_MOTHER_ID": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.PARTICLE_ID)) ) @ F.GET_GENERATION(1),
        "ConstituentMC_GD_MOTHER_Key": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.OBJECT_KEY)) ) @ F.GET_GENERATION(1),
        "ConstituentTRUEP": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.P) ) @ F.GET_GENERATION(1),
        "ConstituentTRUEPT": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PT) ) @ F.GET_GENERATION(1),
        "ConstituentTRUEPX": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PX) ) @ F.GET_GENERATION(1),
        "ConstituentTRUEPY": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PY) ) @ F.GET_GENERATION(1),
        "ConstituentTRUEPZ": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PZ) ) @ F.GET_GENERATION(1),
        "ConstituentTRUEENERGY": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.ENERGY) ) @ F.GET_GENERATION(1),
        "ConstituentQ": F.MAP(F.CHARGE) @ F.GET_GENERATION(1),
        "ConstituentKey": F.MAP(F.OBJECT_KEY) @ F.GET_GENERATION(1),
    })

def JetAdditionalVariables():
    """
    Additional variables for jets.
    """
    return FunctorCollection({
        "Ntracks":      'INFO(9001,-1000)',
        "N90":          'INFO(9002,-1000)',
        "MTF":          'INFO(9003,-1000)',
        "NSatCalo":     'INFO(9004,-1000)',
        "NHasPV":       'INFO(9005,-1000)',
        "CPF":          'INFO(9006,-1000)',
        "JetWidth":     'INFO(9007,-1000)',
        "NSatECAL":     'INFO(9008,-1000)',
        "NSatHCAL":     'INFO(9009,-1000)',
        "NIPChi2Inf4":  'INFO(9010,-1000)',
        "MPT":          'INFO(9011,-1000)',
        "MNF":          'INFO(9012,-1000)',
        "JetWidthNorm": 'INFO(9013,-1000)',
        "vtx_x":        'INFO(9016,-1000)',
        "vtx_y":        'INFO(9017,-1000)',
        "vtx_z":        'INFO(9018,-1000)',
        "JEC":          'INFO(9101,-1000)',
        "NPVsForJEC":   'INFO(9102,-1000)',
        "JECError":     'INFO(9103,-1000)'
    })


def output_tr(**outputs):
    return {"OutputLocation": [v for (k, v) in outputs.items()]}

def make_jet_vars(data_, mctruth_, is_simulation=False):

    jets = find_in_decay(data_, id="CELLjet")


    jet_functors = FunctorCollection({})

    if is_simulation:
        MCjets = make_jets(pt_min = JET_PT_MIN, useMC=True)
        jet_functors += matchjets(jets, MCjets)
        #MCParticles = get_mc_particles('/Event/MC/Particles')
        #pions = find_in_decay(MCParticles, id="pi+")
        #MCParticles = mc_unpackers()["MCParticles"]
        #jet_functors += matchjets(jets, MCParticles)
        #mctruth = MCTruthAndBkgCat(input_particles=jets, name="MCTruthAndBkgCat_functor_jets_{hash}")
        #jet_functors += FC.MCHierarchy(mctruth_alg=mctruth)


    return jet_functors+JetAdditionalVariables()+JetConstituents(data_, mctruth_)


def make_MCTruth_vars(data_, name_=""):
    mctruth = MCTruthAndBkgCat(input_particles=data_, name="MCTruthAndBkgCat_functor_"+name_)
    mc_vars = FC.MCKinematics(mctruth_alg=mctruth)
    mc_vars += FC.MCHierarchy(mctruth_alg=mctruth)
    mc_vars += FC.MCPromptDecay(mctruth_alg=mctruth)
    mc_vars += FC.MCVertexInfo(mctruth_alg=mctruth)
    mc_vars += FC.MCPrimaryVertexInfo(mctruth_alg=mctruth)
    mc_vars += FunctorCollection({"BKGCAT": mctruth.BkgCat})
    return mc_vars, mctruth


def Hlt1Lines_JpsiToMuMu():

    Hlt1_decisions = [
        "Hlt1GlobalDecision", 
        "Hlt1PhysDecision",
        'Hlt1TrackMVADecision',
        'Hlt1TrackMuonMVADecision'
        'Hlt1TwoTrackMVADecision',
        'Hlt1DiMuonHighMassDecision',
        'Hlt1DiMuonLowMassDecision',
        'Hlt1DiMuonSoftDecision',
        "Hlt1DiMuonNoIPDecision",
        'Hlt1SingleHighPtMuonDecision',
        'Hlt1LowPtMuonDecision',
        'Hlt1LowPtDiMuonDecision',
        'Hlt1TrackMuonMVADecision',
        "Hlt1DiMuonNoIP_SSDecision",
        "Hlt1DetJpsiToMuMuPosTagLineDecision",
        "Hlt1DetJpsiToMuMuNegTagLineDecision",
        'Hlt1D2KKDecision',
        'Hlt1D2KPiDecision',
        'Hlt1D2PiPiDecision',
        'Hlt1KsToPiPiDecision',
        ]
    return Hlt1_decisions

def Hlt2Lines_JpsiToMuMu():

    Hlt2_decisions = ["Hlt2_JpsiToMuMuDetachedFullDecision"]

    return Hlt2_decisions

def Hlt2Lines_Wmumuqq():

    Hlt2_decisions = [
        "Hlt2QEE_SingleHighPtMuonFullDecision",
        "Hlt2QEE_SingleHighPtMuonIsoFullDecision",
        "Hlt2QEE_SingleHighPtMuonNoMuIDFullDecision",
        "Hlt2QEE_SingleVHighPtElectronFullDecision",
        "Hlt2QEE_ZToMuMu_DoubleNoMuIDFullDecision",
        "Hlt2QEE_ZToMuMu_SingleNoMuIDFullDecision",
        ]

    return Hlt2_decisions


def SpruceLines_JpsiToMuMu():
    SpruceLines_lines = [
        'SpruceBandQ_TbcToJpsiDz',
        'SpruceBandQ_TbcToJpsiDpKm',
        'SpruceBandQ_TbcToJpsiDzKmPip',
        'SpruceBandQ_TbcToJpsiDzKs',
        'SpruceBandQ_JpsiToMuMuDetached',
        'SpruceBandQ_JpsiToMuMuTightPrompt',
    ]

    return SpruceLines_lines

def make_evt_vars(pvs, Hlt1_decisions_=None, Hlt2_decisions_=None, is_simulation=True, Spruce_decisions_=None):

    odin = get_odin()
    rec_sum=get_rec_summary()

    evt_vars = FunctorCollection({})
    evt_vars += FC.EventInfo()

    if Hlt1_decisions_ is not None:
        evt_vars += FC.SelectionInfo(selection_type="Hlt1", trigger_lines=Hlt1_decisions_)
    if Hlt2_decisions_ is not None:
        evt_vars += FC.SelectionInfo(selection_type="Hlt2", trigger_lines=Hlt2_decisions_)
    if Spruce_decisions_ is not None:
        evt_vars += FC.SelectionInfo(selection_type="Spruce", trigger_lines=Spruce_decisions_)

    if odin:
        evt_vars.update({"EVENTTYPE": F.EVENTTYPE(odin)})
    evt_vars.update({"PV_SIZE": F.SIZE(pvs)})

    evt_vars += FunctorCollection({
        "ALLPVX[nPVs]": F.ALLPVX(pvs),
        "ALLPVY[nPVs]": F.ALLPVY(pvs),
        "ALLPVZ[nPVs]": F.ALLPVZ(pvs),
        "nPVs": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nPVs"),
        "nTTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nTTracks"),
        "nLongTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nLongTracks"),
        "nDownstreamTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nDownstreamTracks"),
        "nUpstreamTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nUpstreamTracks"),
        "nVeloTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nVeloTracks"),
        "nBackTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nBackTracks"),
        "nRich1Hits": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nRich1Hits"),
        "nRich2Hits": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nRich2Hits"),
        "nVPClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nVPClusters"),
        "nUTClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nUTClusters"),
        "nFTClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nFTClusters"),
        "nEcalClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nEcalClusters"),
        "eCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"eCalTot"),
        "hCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"hCalTot"),
    })

    #muon_hits = make_muon_hits()
    
    #evt_vars +=  FunctorCollection({
    #    f"NHITSINMUON_M{station + 1}_R{region}": F.NHITSINMUON(
    #        muon_hits, station=station, region=region
    #   )
    #   for station in [0, 1, 2, 3]
    #   for region in [0, 1, 2, 3]
    #})

    if not is_simulation:
        evt_vars  += FC.LHCInfo()

    return evt_vars

def make_mc_variables(input_data):
    mc_vars =  FunctorCollection({
        "Key": F.OBJECT_KEY,
        'ETA': F.ETA,
        'PHI': F.PHI,
        "THETA": F.THETA,
        'ORIGIN_VX': F.ORIGIN_VX,
        'ORIGIN_VY': F.ORIGIN_VY,
        'ORIGIN_VZ': F.ORIGIN_VZ,
        'MC_LIFETIME': F.MC_LIFETIME,
        'END_VX': F.END_VX,
        'END_VY':F.END_VY,
        'END_VZ':F.END_VZ,
        'PT': F.PT,
        'P': F.P,
        "PX": F.PX,
        "PY": F.PY,
        "PZ": F.PZ,
        "M": F.MASS,
        "ID": F.PARTICLE_ID,
        "Q": F.CHARGE,
        })

    reconstructible_alg_ = MCRectible(input_mctrackinfo=get_mc_track_info())
    mc_vars += FC.MCReconstructible(mcreconstructible_alg=reconstructible_alg_, extra_info=True)
    mcreconstructed_alg_ = MCRected(input_mcparticles=input_data)
    mc_vars += FC.MCReconstructed(mcreconstructed_alg=mcreconstructed_alg_, extra_info=True)

    return mc_vars

def make_Hlt1or2TisTos_vars(Hlt_decisions_, data_, selection_type_ = "Hlt2"):
    return FC.HltTisTos(
        selection_type=selection_type_,
        trigger_lines=Hlt_decisions_,
        data=data_
        )

def make_HltTisTos_vars(data_, Hlt1_decisions_=None, Hlt2_decisions_=None):
    Hlt_vars = FunctorCollection({})

    if Hlt1_decisions_ is not None:
        Hlt_vars += make_Hlt1or2TisTos_vars(Hlt1_decisions_, data_, "Hlt1")
    if Hlt2_decisions_ is not None:
        Hlt_vars += make_Hlt1or2TisTos_vars(Hlt2_decisions_, data_, "Hlt2")

    return Hlt_vars
    
    

