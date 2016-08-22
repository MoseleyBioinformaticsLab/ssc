#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc (Spin System Creator) command-line interface

Usage:
    ssc -h | --help
    ssc --version
    ssc group (--plpath=<path>) (--plformat=<format>) (--stype=<type>) (--dims=<labels>) (--rdims=<labels>) [--results=<path>] [--crspath=<path>]
    ssc group (--descrfile=<path>) [--results=<path>] [--crspath=<path>]

Options:
    -h, --help                   Show this screen.
    --version                    Show version.
    --verbose                    Print what files are processing.
    --plpath=<path>              Path to peak list.
    --plformat=<format>          Peak list format.
    --stype=<type>               Spectrum type.
    --dims=<labels>              Comma-separated dimension labels.
    --rdims=<labels>             Comma-separated root dimension labels.
    --results=<path>             Path where results will be saved.
    --crspath=<path>             Registration algorithm executable path [default: ssc/bin/calculate_registration]
    --descrfile=<path>           Get arguments from description file.
"""

import json
import docopt
import ssc


def main(cmdargs):

    if cmdargs["group"]:
        if cmdargs["--descrfile"]:
            with open(cmdargs["--descrfile"], "r") as infile:
                args = json.load(infile)

            sscreator = ssc.SpinSystemCreator(peaklistpath=args["PeakListPath"],
                                              plformat=args["PeakListFormat"],
                                              spectrumtype=args["SpectrumType"],
                                              dimlabels=args["DimensionLabels"],
                                              rootdims=args["RootDimensions"],
                                              regalgpath=cmdargs["--crspath"],
                                              outputdirpath=cmdargs["--results"])
        else:
            dims = cmdargs["--dims"].split(",") if cmdargs["--dims"] else []
            rdims = cmdargs["--rdims"].split(",") if cmdargs["--rdims"] else []

            sscreator = ssc.SpinSystemCreator(peaklistpath=cmdargs["--plpath"],
                                              plformat=cmdargs["--plformat"],
                                              spectrumtype=cmdargs["--stype"],
                                              dimlabels=dims,
                                              rootdims=rdims,
                                              regalgpath=cmdargs["--crspath"],
                                              outputdirpath=cmdargs["--results"])
        sscreator.group()

args = docopt.docopt(__doc__)
main(args)


# python3 ssc group --plpath=datasets/jr19_hncocacb.pks --plformat=autoassign --stype=HNcoCACB --dims=HN,N,CA/CB-1 --rdims=HN,N
# python3 ssc.pyz group --plpath=datasets/jr19_hncocacb.pks --plformat=autoassign --stype=HNcoCACB --dims=HN,N,CA/CB-1 --rdims=HN,N

# python3 ssc group --descrfile=ProblemDescriptionHNcoCACB_jr19.json