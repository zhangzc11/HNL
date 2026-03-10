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
#from Hlt2Conf.standard_jets import make_particleflow
from Hlt2Conf.standard_jets import make_charged_particles
from RecoConf.standard_particles import make_photons_PF, make_merged_pi0s, get_all_track_selector, standard_protoparticle_filter
from RecoConf.reconstruction_objects import make_pvs
from RecoConf.algorithms_thor import ParticleFilter
from RecoConf.ttrack_selections_reco import make_ttrack_protoparticles
from RecoConf import track_refitting
from PyConf.reading import  reconstruction
from RecoConf.protoparticles import  make_charged_protoparticles

from GaudiKernel.SystemOfUnits import GeV
from Hlt2Conf.standard_jets import build_jets
from PyConf.Algorithms import ParticleFlowMakerMC, ParticleFlowMaker, ParticleMakerForParticleFlow
from RecoConf.data_from_file import mc_unpackers
from PyConf.reading import get_mc_particles

JET_PT_MIN = 10.0 * GeV


def make_comb_mass(nC=1):
    E1 = F.VALUE_OR(0) @ F.CHILD(1, F.VALUE_OR(0) @ F.ENERGY)
    Px1 = F.VALUE_OR(0) @ F.CHILD(1, F.VALUE_OR(0) @ F.PX)
    Py1 = F.VALUE_OR(0) @ F.CHILD(1, F.VALUE_OR(0) @ F.PY)
    Pz1 = F.VALUE_OR(0) @ F.CHILD(1, F.VALUE_OR(0) @ F.PZ)
    E2 = F.VALUE_OR(0) @ F.CHILD(3, F.VALUE_OR(0) @ F.CHILD(nC, F.VALUE_OR(0) @ F.ENERGY))
    Px2 = F.VALUE_OR(0) @ F.CHILD(3, F.VALUE_OR(0) @ F.CHILD(nC, F.VALUE_OR(0) @ F.PX))
    Py2 = F.VALUE_OR(0) @ F.CHILD(3, F.VALUE_OR(0) @ F.CHILD(nC, F.VALUE_OR(0) @ F.PY))
    Pz2 = F.VALUE_OR(0) @ F.CHILD(3, F.VALUE_OR(0) @ F.CHILD(nC, F.VALUE_OR(0) @ F.PZ))
    return  F.VALUE_OR(0) @ fmath.sqrt( (E1+E2)*(E1+E2) - (Px1+Px2)*(Px1+Px2) - (Py1+Py2)*(Py1+Py2) - (Pz1+Pz2) * (Pz1+Pz2) ) 


def make_B_comb_vars():
    return FunctorCollection({
        "MmuBmuN": make_comb_mass(1),
        #"MmuWJet1": make_comb_mass(2),
    })

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
        "OWNPVIP": F.OWNPVIP,           # Impact parameter wrt own PV
        "OWNPVIPCHI2": F.OWNPVIPCHI2,   # Impact parameter chi2 wrt own PV
        "OWNPVX": F.OWNPVX,            # x-coordinate of best PV
        "OWNPVY": F.OWNPVY,            # y-coordinate of best PV
        "OWNPVZ": F.OWNPVZ,
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
        "VTXCHI2": F.CHI2,         # Vertex fit chi2/ndf
        "VTXCHI2NDOF": F.CHI2DOF,         # Vertex fit chi2/ndf
        "END_VX": F.END_VX,               # x-coordinate of decay vertex
        "END_VY": F.END_VY,               # y-coordinate of decay vertex
        "END_VZ": F.END_VZ,               # z-coordinate of decay vertex
        #"END_VR": fmath.sqrt(F.END_VZ*F.END_VZ + F.END_VY*F.END_VY + F.END_VX*F.END_VX),
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
        "SDOCA13": F.SDOCA(1, 3),
        "DOCA13CHI2": F.DOCACHI2(1, 3),
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
        "ISMUONTIGHT": F.ISMUONTIGHT,              # Boolean: is it identified as a muon 0 or 1?
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

def make_ttrack():
    track_fitter=track_refitting.TRACK_FIT_TYPE_PRKALMAN
    raw_ttrack = reconstruction()["Ttracks"]
    fitted_ttrack = track_refitting.refit_tracks(
            raw_ttrack,
            get_clusters_from_track=True,
            track_fitter=track_fitter,
            update_ghost_prob=True,
            disable_UT=False,
            )
    return fitted_ttrack


def make_ttrack_protoparticles_diy():
    return make_charged_protoparticles(tracks = make_ttrack(),
                                       rich_pids = None,
                                       calo_pids = None,
                                       muon_pids = None,
                                       track_types = ["Ttrack"]
                                       )

def make_charged_particles_Ttrack(
    get_track_selector=get_all_track_selector,
    make_protoparticle_filter=standard_protoparticle_filter,
    c_over_e_cut=0,
    ):
    chargedProtos = make_ttrack_protoparticles_diy()
    return ParticleMakerForParticleFlow(
        InputProtoParticles=chargedProtos,
        TrackPredicate=get_track_selector(),
        c_over_e_cut=c_over_e_cut,
        ProtoParticlePredicate=make_protoparticle_filter(),
        PrimaryVertices=make_pvs(),
    ).Output


def make_particleflow(name='PF_{hash}'):
    charged_long = ParticleFilter(
        make_charged_particles(track_type='Long'),
        F.FILTER(F.require_all( F.ISMUON != 1, F.IS_ABS_ID("mu-") != 1, fmath.abs(F.PARTICLE_ID) != 13))
        )

    #charged_down = make_charged_particles(track_type='Downstream')

    charged_down = ParticleFilter(
        make_charged_particles(track_type='Downstream'),
        F.FILTER(F.require_all( F.ISMUON != 1, F.IS_ABS_ID("mu-") != 1, fmath.abs(F.PARTICLE_ID) != 13))
        )

    #charged_T = make_charged_particles_Ttrack()

    photons = make_photons_PF()
    pi0s = make_merged_pi0s()

    return ParticleFlowMaker(
        #Inputs=[charged_long, charged_down, charged_T, photons, pi0s],
        Inputs=[charged_long, charged_down, photons, pi0s],
        name=name).Output

def make_particleflowMC():
    #requestedParticlesPIDs=[11, -11, 13, -13, 211, -211, 321, -321, 2112, -2112, 2212, -2212, 22, 310, 311, 130, 111, 3122, -3122, 433, -433, 20433, -20433, 20433], 
    requestedParticlesPIDs=[11, -11, 13, -13, 22, 111, 113, 115, 130, -211, 211, -213, 213, 215, 221, 223, 225, 310, -311, 311, -313, 313, 315, -321, 321, -323, 323, 331, 333, -411, -413, 413, -415, -421, -423, -425, -431, -433, 1114, -2112, 2112, 2114, -2212, 2212, -2214, 2214, 2224, 3112, -3122, 3122, -3212, 3214, -3222, 3224, -4122, -4222, 10113, -10211, -10213, 10213, 10223, 10311, 10313, 10321, -10413, -10433, -20213, 20213, 20223, 20323, -20413, -20433]
    banPIDs=[12, -12, 14, -14, 16, -16, 23] # ν_e, ν̄_e, ν_μ, ν̄_μ, ν_τ, ν̄_τ, Z0
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
        Cut=F.require_all(F.DR2 < (1.0**2)),
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
        "MCClosestJet_PX":  F.VALUE_OR(-99) @ F.PX @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_PY":  F.VALUE_OR(-99) @ F.PY @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_PZ":  F.VALUE_OR(-99) @ F.PZ @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_ETA":  F.VALUE_OR(-99) @ F.ETA @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_PHI":  F.VALUE_OR(-99) @ F.PHI @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_E":  F.VALUE_OR(-99) @ F.ENERGY @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        "MCClosestJet_ID":  F.VALUE_OR(-99) @ F.PARTICLE_ID @ F.TO @ F.ENTRY_WITH_MIN_REL_VALUE_OF( DR2 ).bind( F.RELATIONS.bind(F.TES(alg.OutputRelations), F.FORWARDARGS), F.FORWARDARGS ),
        }
    )

    return functors

def JetConstituents(data__, MCTRUTH):
    """
    Get the constituents of the jets.
    """
    var_const = FunctorCollection({
        "nCons": F.NINGENERATION(F.ALL,1),
        "nChargedCons": F.NINGENERATION(F.CHARGE != 0, 1),
        "nNeutralCons": F.NINGENERATION(F.CHARGE == 0, 1),
        # Store an array of jet constituents' Px, Py, Pz
        "ConPX": F.MAP(F.PX) @ F.GET_GENERATION(1),
        "ConPY": F.MAP(F.PY) @ F.GET_GENERATION(1),
        "ConPZ": F.MAP(F.PZ) @ F.GET_GENERATION(1),
        "ConPT": F.MAP(F.PT) @ F.GET_GENERATION(1),
        "ConETA": F.MAP(F.ETA) @ F.GET_GENERATION(1),
        "ConPHI": F.MAP(F.PHI) @ F.GET_GENERATION(1),
        "ConE": F.MAP(F.ENERGY) @ F.GET_GENERATION(1),
        "ConMass": F.MAP(F.MASS) @ F.GET_GENERATION(1),
        "ConPID": F.MAP(F.PARTICLE_ID) @ F.GET_GENERATION(1),
        "ConNHITS": F.MAP( F.VALUE_OR(0) @ F.NHITS @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumNHITS":  F.SUM (F.VALUE_OR(0) @ F.NHITS @ F.TRACK),
        "ConNFTHITS": F.MAP( F.VALUE_OR(0) @ F.NFTHITS @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNFTHITS": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NFTHITS @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNVPHITS": F.MAP( F.VALUE_OR(0) @ F.NVPHITS @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNVPHITS": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NVPHITS @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNVPHITSC": F.MAP( F.VALUE_OR(0) @ F.NVPHITSC @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNVPHITSC": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NVPHITSC @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNVPHITSA": F.MAP( F.VALUE_OR(0) @ F.NVPHITSA @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNVPHITSA": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NVPHITSA @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNVPHOLES": F.MAP( F.VALUE_OR(0) @ F.NVPHOLES @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNVPHOLES": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NVPHOLES @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNFTHOLES": F.MAP( F.VALUE_OR(0) @ F.NFTHOLES @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNFTHOLES": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NFTHOLES @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNVPLAYERS": F.MAP( F.VALUE_OR(0) @ F.NVPLAYERS @ F.TRACK ) @ F.GET_GENERATION(1),
        #"ConSumNVPLAYERS": F.SUM @ F.MAP( F.VALUE_OR(0) @ F.NVPLAYERS @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConNUTHITS": F.MAP( F.VALUE_OR(0) @ F.NUTHITS @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumNUTHITS": F.COUNT_IF( ( F.VALUE_OR(0) @ F.NUTHITS @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKHASVELO": F.MAP( F.VALUE_OR(0) @ F.TRACKHASVELO @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKHASVELO": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKHASVELO @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKHASUT": F.MAP( F.VALUE_OR(0) @ F.TRACKHASUT @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKHASUT": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKHASUT @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKHAST": F.MAP( F.VALUE_OR(0) @ F.TRACKHAST @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKHAST": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKHAST @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKISLONG": F.MAP( F.VALUE_OR(0) @ F.TRACKISLONG @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKISLONG": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKISLONG @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKISDOWNSTREAM": F.MAP( F.VALUE_OR(0) @ F.TRACKISDOWNSTREAM @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKISDOWNSTREAM": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKISDOWNSTREAM @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKISTTRACK": F.MAP( F.VALUE_OR(0) @ F.TRACKISTTRACK @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKISTTRACK": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKISTTRACK @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKISUPSTREAM": F.MAP( F.VALUE_OR(0) @ F.TRACKISUPSTREAM @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKISUPSTREAM": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKISUPSTREAM @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),
        "ConTRACKISVELO": F.MAP( F.VALUE_OR(0) @ F.TRACKISVELO @ F.TRACK ) @ F.GET_GENERATION(1),
        "ConSumTRACKISVELO": F.COUNT_IF( ( F.VALUE_OR(0) @ F.TRACKISVELO @ F.TRACK) == 1 ) @ F.GET_GENERATION(1),

        "ConINECAL": F.MAP( F.VALUE_OR(0) @ F.INECAL ) @ F.GET_GENERATION(1),
        "ConSumINECAL": F.COUNT_IF( (F.VALUE_OR(0) @ F.INECAL) == 1 ) @ F.GET_GENERATION(1),
        "ConINHCAL": F.MAP( F.VALUE_OR(0) @ F.INHCAL ) @ F.GET_GENERATION(1),
        "ConSumINHCAL": F.COUNT_IF( (F.VALUE_OR(0) @ F.INHCAL) == 1 ) @ F.GET_GENERATION(1),
        "ConINBREM": F.MAP( F.VALUE_OR(0) @ F.INBREM ) @ F.GET_GENERATION(1),
        "ConSumINBREM": F.COUNT_IF( (F.VALUE_OR(0) @ F.INBREM) == 1 ) @ F.GET_GENERATION(1),
        "ConINMUON": F.MAP( F.VALUE_OR(0) @ F.INMUON ) @ F.GET_GENERATION(1),
        "ConSumINMUON": F.COUNT_IF( (F.VALUE_OR(0) @ F.INMUON) == 1 ) @ F.GET_GENERATION(1),
        #"ConISMUON": F.MAP( F.VALUE_OR(0) @ F.ISMUON ) @ F.GET_GENERATION(1), # doesn't work for 2J??
        #"ConSumISMUON": F.COUNT_IF( (F.VALUE_OR(0) @ F.ISMUON) == 1 ) @ F.GET_GENERATION(1),
        "ConQ": F.MAP(F.CHARGE) @ F.GET_GENERATION(1),
        "ConKey": F.MAP(F.OBJECT_KEY) @ F.GET_GENERATION(1),
    })
    if MCTRUTH is not None:
        var_const += FunctorCollection({
            "ConTRUEID": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.PARTICLE_ID) ) @ F.GET_GENERATION(1),
            "ConTRUEID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.PARTICLE_ID)) == 25) @ F.GET_GENERATION(1),
            "ConMC_MOTHER_ID": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(1, F.PARTICLE_ID)) ) @ F.GET_GENERATION(1),
            "ConMC_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(1, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_MOTHER_Key": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(1, F.OBJECT_KEY)) ) @ F.GET_GENERATION(1),
            "ConMC_GD_MOTHER_ID": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.PARTICLE_ID)) ) @ F.GET_GENERATION(1),
            "ConMC_GD_MOTHER_ID_N511": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.PARTICLE_ID))) == 511) @ F.GET_GENERATION(1),
            "ConMC_GD2_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD3_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(3, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD4_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(4, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD5_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(5, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD6_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(6, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD7_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(7, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD8_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(8, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD9_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(9, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GD10_MOTHER_ID_N25": F.COUNT_IF( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(10, F.PARTICLE_ID))) == 25) @ F.GET_GENERATION(1),
            "ConMC_GDall_MOTHER_ID_N25": F.COUNT_IF( 
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(1, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(3, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(4, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(5, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(6, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(7, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(8, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(9, F.PARTICLE_ID))) == 25 ) |
            ( fmath.abs(F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(10, F.PARTICLE_ID))) == 25 )
            ) @ F.GET_GENERATION(1),
            "ConMC_GD_MOTHER_Key": F.MAP( F.VALUE_OR(0) @ MCTRUTH(F.MC_MOTHER(2, F.OBJECT_KEY)) ) @ F.GET_GENERATION(1),
            "ConTRUEP": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.P) ) @ F.GET_GENERATION(1),
            "ConTRUEPT": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PT) ) @ F.GET_GENERATION(1),
            "ConTRUEPX": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PX) ) @ F.GET_GENERATION(1),
            "ConTRUEPY": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PY) ) @ F.GET_GENERATION(1),
            "ConTRUEPZ": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.PZ) ) @ F.GET_GENERATION(1),
            "ConTRUEE": F.MAP( F.VALUE_OR(-1) @ MCTRUTH(F.ENERGY) ) @ F.GET_GENERATION(1),
    })

    return var_const

def JetAdditionalVariables():
    """
    Additional variables for jets.
    """
    return FunctorCollection({
        "Ntracks":      'INFO(9001,-1000)', # number of charged components
        "N90":          'INFO(9002,-1000)', # 
        "MTF":          'INFO(9003,-1000)', # max pt of charged components / jet pt
        "NSatCalo":     'INFO(9004,-1000)',
        "NHasPV":       'INFO(9005,-1000)',
        "CPF":          'INFO(9006,-1000)', # vector sum of charged components' pt / jet pt
        "JetWidth":     'INFO(9007,-1000)', # DR between components and jet * component's pt
        "NSatECAL":     'INFO(9008,-1000)', 
        "NSatHCAL":     'INFO(9009,-1000)',
        "NIPChi2Inf4":  'INFO(9010,-1000)',
        "MPT":          'INFO(9011,-1000)', # max pt of charged components
        "MNF":          'INFO(9012,-1000)', # max pt of neutral componnets / jet pt
        "JetWidthNorm": 'INFO(9013,-1000)', # JetWidth/sum of component's pt
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


def Hlt1Lines_Bmumuqq():

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
        ]
    return Hlt1_decisions

def Hlt2Lines_Bmumuqq():

    Hlt2_decisions = [
        "Hlt2QEE_SingleHighPtMuonFullDecision",
        "Hlt2QEE_SingleHighPtMuonIsoFullDecision",
        "Hlt2QEE_SingleHighPtMuonNoMuIDFullDecision",
        "Hlt2QEE_SingleVHighPtElectronFullDecision", # 35 GeV
        "Hlt2QEE_SingleHighPtElectronFullDecision", # 17.5 GeV
        "Hlt2QEE_ZToMuMu_DoubleNoMuIDFullDecision",
        "Hlt2QEE_ZToMuMu_SingleNoMuIDFullDecision",
        "Hlt2QEE_IncJet15GeVFullDecision",
        "Hlt2QEE_IncJet25GeVFullDecision",
        "Hlt2QEE_IncJet35GeVFullDecision",
        "Hlt2QEE_IncJet45GeVFullDecision",
        "Hlt2QEE_DiSVTagJet10GeVFullDecision",
        "Hlt2QEE_DiSVTagJet15GeVFullDecision",
        "Hlt2QEE_DiSVTagJet20GeVFullDecision",
        "Hlt2QEE_DiSVTagJet25GeVFullDecision",
        "Hlt2QEE_DiSVTagJet30GeVFullDecision",
        "Hlt2QEE_DiSVTagJet35GeVFullDecision",
        "Hlt2QEE_DiTopoTagJet10GeVFullDecision",
        "Hlt2QEE_DiTopoTagJet15GeVFullDecision",
        "Hlt2QEE_DiTopoTagJet20GeVFullDecision",
        "Hlt2QEE_DiTopoTagJet25GeVFullDecision",
        "Hlt2QEE_DiTopoTagJet30GeVFullDecision",
        "Hlt2QEE_DiTopoTagJet35GeVFullDecision",
        "Hlt2QEE_DiJetIncSVTag_pT25M40FullDecision",
        ]

    return Hlt2_decisions


def SpruceLines_Bmumuqq():
    SpruceLines_lines = [
        'SpruceQEE_SingleHighPtMuon',  # 15 GeV, no prescale
        'SpruceQEE_SingleHighPtMuonIso', #12.5 GeV, prescale = 1.0*0.1, Hlt2*Spruce
        'SpruceQEE_SingleJet15', # prescale = 0.05 * 0.05
        'SpruceQEE_SingleJet25', # prescale = 0.1 * 0.05
        'SpruceQEE_SingleJet35', # prescale = 0.5 * 0.05
        'SpruceQEE_SingleJet45', # prescale = 1.0 * 0.05
        'SpruceQEE_diSVTag1010', # prescale = 0.0025 * 1.0
        'SpruceQEE_diSVTag1515', # prescale = 0.01 * 1.0
        'SpruceQEE_diSVTag2020', # prescale = 0.05 * 1.0
        'SpruceQEE_diSVTag2525', # prescale = 0.1 * 1.0
        'SpruceQEE_diSVTag3030', # prescale = 0.25 * 1.0
        'SpruceQEE_diSVTag3535', # prescale = 1.0 * 1.0
        'SpruceQEE_diTopoTag1010', # prescale = 0.001 * 1.0
        'SpruceQEE_diTopoTag1515', # prescale = 0.001 * 1.0
        'SpruceQEE_diTopoTag2020', # prescale = 0.005 * 1.0
        'SpruceQEE_diTopoTag2525', # prescale = 0.05 * 1.0
        'SpruceQEE_diTopoTag3030', # prescale = 0.05 * 1.0
        'SpruceQEE_diTopoTag3535', # prescale = 0.1 * 1.0
        'SpruceQEE_Dijets1515', # prescale = 0.1 * 0.05
        'SpruceQEE_Dijets2020', # prescale = 0.25 * 0.05
        'SpruceQEE_Dijets2525', # prescale = 0.5 * 0.05
        'SpruceQEE_Dijets3030', # prescale = 0.75 * 0.05
        'SpruceQEE_Dijets3535', # prescale = 1.0 * 0.05
        'SpruceQEE_IncSVTagDijets', # no prescale
        'SpruceQEE_WJet', # muon pt > 20GeV, jet pt > 10 GeV
        'SpruceQEE_WJetJet', # muon pt > 20GeV, jet pt > 10 GeV
        'SpruceQEE_SingleHighPtElectron', # no prescale
        'SpruceQEE_SingleVHighPtElectron', # prescale = 1.0 *0.05
        'SpruceQEE_ZToMuMu', # no prescale
        'SpruceQEE_ZToMuMu_SingleNoMuID', # no prescale
        'SpruceQEE_ZToMuMu_DoubleNoMuID', # no prescale
        'SpruceBandQ_DiMuonInc',
        'SpruceBandQ_DiMuonSameSignInc',
        'SpruceBandQ_DiMuonSameSignHighMass',
        'SpruceBandQ_DiMuonIncHighPT',
        'SpruceBandQ_DiMuonSameSignIncHighPT',
        'SpruceBandQ_DiMuonSoft',
        'SpruceBandQ_BuForSpectroscopy',
        'SpruceBandQ_BdForSpectroscopy',
        'SpruceBandQ_BsForSpectroscopy',
        'SpruceBandQ_BcForSpectroscopy',
        'SpruceBandQ_BudForSpectroscopySL',
        'SpruceBandQ_BsForSpectroscopySL',
        'SpruceBandQ_BcForSpectroscopySL',
        'SpruceBandQ_JpsiToMuMuDetached',
        'SpruceBandQ_JpsiToMuMuTightPrompt',
        'SpruceBandQ_Psi2SToMuMuDetached',
        'SpruceBandQ_Psi2SToMuMuTightPrompt',
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
        #"eCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"eCalTot"),
        #"hCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"hCalTot"), # doesn't work for qee data
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
        "END_VR": fmath.sqrt(F.END_VZ*F.END_VZ + F.END_VY*F.END_VY + F.END_VX*F.END_VX),
        'PT': F.PT,
        'P': F.P,
        "PX": F.PX,
        "PY": F.PY,
        "PZ": F.PZ,
        "M": F.MASS,
        "E": F.ENERGY,
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
    
    
