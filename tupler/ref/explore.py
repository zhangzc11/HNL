###############################################################################
# (c) Copyright 2025 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################

# Imports
import sys

import GaudiPython as GP
from Moore import options
from PyConf.application import configure, configure_input
from PyConf.reading import tes_root_mc

LHCb = GP.gbl.LHCb
import argparse


def printMCParticles(evt):
    n_size = evt['/Event/HLT2/MC/Particles'].size()
    for ip in range(n_size):
        print(evt['/Event/HLT2/MC/Particles'][ip])


def advance_decision(decision, source="Hlt2"):
    """
    Advance to the next event with a positive decision for the given HLT
    line.

    Args:
        decision (str): Name of the HLT line to select events on.
        source (str, optional): Selection stage. Defaults to "Hlt2".

    Raises:
        SystemExit: If an event with a positive decision is not found.
    """

    loc = f"/Event/{source}/DecReports"
    while True:
        appMgr.run(1)
        if not evt["/Event"]:
            sys.exit("Did not find positive {0} decision".format(decision))
        reports = evt[loc]
        report = reports.decReport("{0}Decision".format(decision))
        if report.decision() == 1:
            break


def list_fired_triggers(source="Hlt2"):
    """
    Print the names of the HLT lines that have a positive decision.

    Args:
        source (str, optional): Selection stage. Defaults to "Hlt2".
    """

    loc = f"/Event/{source}/DecReports"
    reports = evt[loc]
    for i in reports.decisionNames():
        if reports.decReport(i).decision() == 1:
            print(i)


# Argument parser
parser = argparse.ArgumentParser(
    usage="./run python -i %(prog)s --input_file xxx.[dst, mdf] (-i process) (-s stream) (--sim)",
    description="Inspect Moore output",
)
parser.add_argument("-i", "--input_file", type=str, help="Input file", required=True)
parser.add_argument(
    "-p",
    "--input_process",
    type=str,
    help="Last process run - 'Turbo' or 'Spruce' or 'Hlt2' or 'TurboSpruce' or `TurboPass`",
    required=False,
    default="Turbo",
)
parser.add_argument(
    "-s",
    "--input_stream",
    type=str,
    help="Stream if not Hlt2 input_process",
    required=False,
    default="DAQ",
)
parser.add_argument(
    "--sim",
    action="store_true",
    help="Run in simulation mode (default: False)",
    # 'default=False' is implied by store_true
)

args = parser.parse_args()

options.data_type = "Upgrade"
options.simulation = args.sim
options.input_files = [args.input_file]
options.input_type = "ROOT" if "dst" in args.input_file else "RAW"
options.root_ioalg_name = "RootIOAlgExt"
options.evt_max = -1
options.gaudipython_mode = True
options.input_stream = args.input_stream
options.geometry_version = "run3/trunk"
options.conditions_version = "master"
options.dstformat = "uDST"
options.input_process = args.input_process

config = configure_input(options)

# Following is needed as binds from https://gitlab.cern.ch/lhcb/Moore/-/blob/master/Hlt/Moore/python/Moore/LbExec.py#L41
# are not included when using GP
tes_root_mc.global_bind(dstformat="uDST", input_process=options.input_process)

from GaudiConf.reading import do_unpacking

algs = do_unpacking(
    input_process=options.input_process,
    has_mc_data=options.simulation,
    include_hlt1=True,
)

config.update(configure(options, algs))
appMgr = GP.AppMgr()
evt = appMgr.evtsvc()
appMgr.run(1)
evt.dump()

#printMCParticles(evt)
