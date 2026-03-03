from FunTuple import FunctorCollection
import FunTuple.functorcollections as FC
import Functors as F
from DecayTreeFitter import DecayTreeFitter
from DaVinciMCTools import MCTruthAndBkgCat
from PyConf.reading import get_rec_summary, get_odin


def make_base_vars(pvs):
    return FunctorCollection({
        "M": F.MASS,
        "P": F.P,
        "E": F.ENERGY,
        "ETA": F.ETA,
        "PHI": F.PHI,
        "PT": F.PT,
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
        "OWNPV_IPCHI2": F.OWNPVIPCHI2,    # Impact parameter Chi2 with respect to the own PV
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
        "MIN_PT": F.MIN(F.PT),            # Minimum PT of daughters
        "SUM_PT": F.SUM(F.PT),            # Sum of daughters' PT
        "MAX_P": F.MAX(F.P),              # Maximum momentum of daughters
        "MIN_P": F.MIN(F.P),              # Minimum momentum of daughters
        "SUM_P": F.SUM(F.P),              # Sum of daughters' momentum
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
        #"TRACKTYPE": F.TRACKTYPE,
        "PIDp": F.PID_P,              # Proton PID likelihood
        "PIDK": F.PID_K,              # Kaon PID likelihood
        "PIDPi": F.PID_PI,            # Pion PID likelihood
        "PIDe": F.PID_E,              # Electron PID likelihood
        "PIDmu": F.PID_MU,            # Muon PID likelihood
        # PROBNNs
        "PROBNN_pi": F.PROBNN_PI,        # Neural net probability of being a pion
        "PROBNN_p": F.PROBNN_P,          # Neural net probability of being a proton
        "PROBNN_K": F.PROBNN_K,          # Neural net probability of being a kaon
        "PROBNN_e": F.PROBNN_E,          # Neural net probability of being an electron
        "PROBNN_mu": F.PROBNN_MU,        # Neural net probability of being a muon
        "PROBNN_GHOST": F.PROBNN_GHOST,  # Neural net probability of being a ghost track
        "GHOSTPROB": F.GHOSTPROB,
        "ISMUON": F.ISMUON,              # Boolean: is it identified as a muon 0 or 1?
        "INMUON": F.INMUON,
        "INECAL": F.INECAL,
        "INHCAL": F.INHCAL,
        "OWNPVIPCHI2": F.OWNPVIPCHI2,
        "MINIPCHI2": F.MINIPCHI2(pvs),
        "HASBREM": F.HASBREM,
        "TRCHI2":    F.CHI2 @ F.TRACK,
        "TRCHI2DOF": F.CHI2DOF @ F.TRACK,
        # Additional track related info
        # F.TRACK gets the track object
        # F.NDOF gets the NDOF for that track
        # F.VALUE_OR(-1) means if no value exists return -1 instead of failing
        "QOVERP": F.QOVERP @ F.TRACK,
        "TRNDOF": F.VALUE_OR(-1) @ F.NDOF @ F.TRACK,                # NDOF in track fit, if this is higher then more hits used in fit
        "NHITS": F.VALUE_OR(-1) @ F.NHITS @ F.TRACK,                # Total number of hits in all detectors
        "NVPHITS": F.VALUE_OR(-1) @ F.NVPHITS @ F.TRACK,            # Total number of hits in VELO phi sensors
        "NVPHITSC":  F.VALUE_OR(-1) @ F.NVPHITSC @ F.TRACK,
        "NUTHITS": F.VALUE_OR(-1) @ F.NUTHITS @ F.TRACK,            # Total number of hits in UT
        "NVPHITSA":  F.VALUE_OR(-1) @ F.NVPHITSA @ F.TRACK,
        "NFTHITS": F.VALUE_OR(-1) @ F.NFTHITS @ F.TRACK,            # Total number of hits in Fibre Tracker (SciFi)
        "TRACKHISTORY": F.VALUE_OR(-1) @ F.TRACKHISTORY @ F.TRACK,  # Track reconstruction history
        "TRACKHASVELO": F.VALUE_OR(-1) @ F.TRACKHASVELO @ F.TRACK,
        "TRACKHASUT": F.VALUE_OR(-1) @ F.TRACKHASUT @ F.TRACK,
        "PPHASRICH": F.VALUE_OR(-1) @ F.PPHASRICH @ F.PROTOPARTICLE,
        "ISLONG": F.VALUE_OR(-1) @ F.TRACKISLONG @ F.TRACK,
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

def make_MCTruth_vars(data_, name_=""):
    mctruth = MCTruthAndBkgCat(input_particles=data_, name="MCTruthAndBkgCat_functor_"+name_)
    mc_vars = FC.MCKinematics(mctruth_alg=mctruth)
    mc_vars += FC.MCHierarchy(mctruth_alg=mctruth)
    mc_vars += FC.MCPrimaryVertexInfo(mctruth_alg=mctruth)
    mc_vars += FunctorCollection({"BKGCAT": mctruth.BkgCat})
    return mc_vars

def make_evt_vars(pvs, line, selection_type="Hlt2", is_simulation=True):

    odin = get_odin()
    rec_sum=get_rec_summary()

    evt_vars = FunctorCollection({})
    evt_vars += FC.EventInfo()

    evt_vars += FC.SelectionInfo(selection_type=selection_type, trigger_lines=[line + "Decision"])

    if odin:
        evt_vars.update({"EVENTTYPE": F.EVENTTYPE(odin)})
    evt_vars.update({"PV_SIZE": F.SIZE(pvs)})

    evt_vars += FunctorCollection({
        "ALLPVX[nPVs]": F.ALLPVX(pvs),
        "ALLPVY[nPVs]": F.ALLPVY(pvs),
        "ALLPVZ[nPVs]": F.ALLPVZ(pvs),
        "numberPVs": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nPVs"),
        "nTTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nTTracks"),
        "nLongTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nLongTracks"),
        "nDownstreamTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nDownstreamTracks"),
        "nUpstreamTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nUpstreamTracks"),
        "nVeloTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nVeloTracks"),
        "nBackTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nBackTracks"),
        "nRich1Hits": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nRich1Hits"),
        "nRich2Hits": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nRich2Hits"),
        "nVPClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nVPClusters"),
        "nFTClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nFTClusters"),
        "eCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"eCalTot"),
        "hCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"hCalTot"),
        "nEcalClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nEcalClusters"),
    })

    evt_vars  += FC.EventInfo()

    if not is_simulation:
        evt_vars  += FC.LHCInfo()

    return evt_vars

def make_mc_variables():
    return FunctorCollection({
        'ETA': F.ETA,
        'PHI': F.PHI,
        'ORIGIN_VX': F.ORIGIN_VX,
        'ORIGIN_VY': F.ORIGIN_VY,
        'ORIGIN_VZ': F.ORIGIN_VZ,
        'END_VX': F.END_VX,
        'END_VY':F.END_VY,
        'END_VZ':F.END_VZ,
        'PT': F.PT,
        'P': F.FOURMOMENTUM,
        "PX": F.PX,
        "PY": F.PY,
        "PZ": F.PZ,
        "M": F.MASS,
        "ID": F.PARTICLE_ID,
        "Q": F.CHARGE,
        })

def make_HltTisTos_vars(Hlt2_decisions_, data_, selection_type_ = "Hlt2"):
    return FC.HltTisTos(
        selection_type=selection_type_,
        trigger_lines=[f"{x}Decision" for x in Hlt2_decisions_],
        data=data_
        )

def make_Hlt2TisTos_vars(data_, Hlt2_decisions_=None):
    Hlt2_decisions = ["Hlt2_JpsiToMuMuDetachedFull"]
    if Hlt2_decisions_ is not None:
        Hlt2_decisions = Hlt2_decisions_
    return make_HltTisTos_vars(Hlt2_decisions, data_, "Hlt2")


