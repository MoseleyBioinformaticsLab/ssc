#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc (Spin System Creator) command-line interface

Usage:
    ssc -h | --help
    ssc --version
    ssc group (--plpath=<path>) (--plformat=<format>) (--stype=<type>) (--dims=<labels>) (--rdims=<labels>) [--result=<path>] [--crs=<path>] [--view]
    ssc visualize <grouping_result> <x_idx> <y_idx> <x_label> <y_label> <plot_title>

Options:
    -h, --help                   Show this screen.
    --version                    Show version.
    --view                       Print results of the grouping.
    --plpath=<path>              Path to peak list.
    --plformat=<format>          Peak list format.
    --stype=<type>               Spectrum type.
    --dims=<labels>              Comma-separated dimension labels.
    --rdims=<labels>             Comma-separated root dimension labels.
    --crs=<path>                 Registration algorithm executable path [default: ssc/bin/CRS_EXE]
    --result=<path>              Path to directory where results will be saved.
"""

from . import docopt
from . import ssc
from . import __version__
from . import visualize


def main(cmdargs):
    if cmdargs["group"]:
        dims = cmdargs["--dims"].split(",") if cmdargs["--dims"] else []
        rdims = cmdargs["--rdims"].split(",") if cmdargs["--rdims"] else []
        view = cmdargs["--view"]

        sscreator = ssc.SpinSystemCreator(peaklist_path=cmdargs["--plpath"],
                                          plformat=cmdargs["--plformat"],
                                          spectrum_type=cmdargs["--stype"],
                                          labels=dims,
                                          root_dims=rdims,
                                          regalg_path=cmdargs["--crs"],
                                          grouping_result_path=cmdargs["--result"])
        sscreator.group(view=view)

    if cmdargs["visualize"]:

        visualize.visualize_clusters(clusters_path=cmdargs["<grouping_result>"],
                                     x_idx=int(cmdargs["<x_idx>"]),
                                     y_idx=int(cmdargs["<y_idx>"]),
                                     x_label=cmdargs["<x_label>"],
                                     y_label=cmdargs["<y_label>"],
                                     title=cmdargs["<plot_title>"])

args = docopt.docopt(__doc__, version=__version__)
main(args)
