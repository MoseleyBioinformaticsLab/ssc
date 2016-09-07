ssc
===

Spin System Creator (ssc)

This package provides `SpinSystemCreator` class that groups peaks from
a single peak list into clusters (groups of peaks that belong to the same
spin system).

.. code:: bash

   ssc (Spin System Creator) command-line interface

   Usage:
       ssc -h | --help
       ssc --version
       ssc group (--plpath=<path>) (--plformat=<format>) (--stype=<type>) (--dims=<labels>) (--rdims=<labels>) [--results=<path>] [--crspath=<path>]
       ssc group (--descrfile=<path>) [--results=<path>] [--crspath=<path>]
       ssc group (--stdin) [--results=<path>] [--crspath=<path>]

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
       --stdin                      Get arguments from standard input.
