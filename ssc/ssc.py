#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.ssc
~~~~~~~

This module provides :class:`~ssc.ssc.SpinSystemCreator` class that groups peaks from
a single peak list into spin systems.
"""

import sys
import os
import datetime
import json

import peaklistparsers as plp
import physicalentities as pe
import registration
import grouping


class SpinSystemCreator(object):
    """Perform internal peak list registration and internal grouping of a single peak list,
    i.e. group peaks that belong to the same spin system into clusters."""

    formats = {"sparky": "-sparky",
               "autoassign": "-autoassign",
               "json": "-json"}

    def __init__(self, peaklistpath, plformat, spectrumtype, dimlabels, rootdims, regalgpath, outputdirpath=None):
        """Initialize SpinSystemCreator.

        :param str peaklistpath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension labels.
        :param list rootdims: List of root dimension labels.
        :param str plformat: Peak list format.
        :param str regalgpath: Path to registration algorithm executable.
        :param str outputdirpath: Path to directory where save the results.
        """
        self.peaklistpath = os.path.normpath(peaklistpath)
        self.spectrumtype = spectrumtype
        self.dimlabels = dimlabels
        self.plformat = plformat
        self.rootdims = rootdims
        self.regalgpath = os.path.normpath(regalgpath)

        if not outputdirpath:
            self.outputdirpath = os.path.join(os.getcwd(),
                                              "results",
                                              datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        else:
            self.outputdirpath = os.path.join(os.path.normpath(outputdirpath),
                                              "results",
                                              datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        self.reg_result_dir = os.path.join(self.outputdirpath, "registration_result")
        self.group_result_dir = os.path.join(self.outputdirpath, "grouping_result")
        self.pl_dir = os.path.join(self.outputdirpath, "peaklists")
        self._makedirs(self.outputdirpath, self.reg_result_dir, self.group_result_dir, self.pl_dir)

    @staticmethod
    def _makedirs(*dirpath):
        """Make directories using directory path.

        :param str dirpath: Path to where make a directory.
        :return: None
        :rtype: None
        """
        for directory in dirpath:
            if not os.path.exists(directory):
                os.makedirs(directory)

    @staticmethod
    def parse(peaklistpath, spectrumtype, dimlabels, plformat):
        """Parse single peak list.

        :param str peaklistpath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension labels.
        :param str plformat: Peak list format.
        :return: Peak list.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        if plformat == "sparky":
            peaklist = plp.SparkyPeakListParser.parse(peaklistpath, spectrumtype, dimlabels, plformat)
        elif plformat == "autoassign":
            peaklist = plp.AutoAssignPeakListParser.parse(peaklistpath, spectrumtype, dimlabels, plformat)
        elif plformat == "json":
            pass
        else:
            raise TypeError('Unknown peak list format: "{}"'.format(plformat))

        return peaklist

    @staticmethod
    def calculate_registration(regalgpath, inputlistpath, rootlistpath, plformat, regresultpath, rootdims, inputdims):
        """Execute registration algorithm binary.

        :param str inputlistpath: Path to input peak list.
        :param str rootlistpath: Path to root peak list.
        :param str plformat: Peak list format.
        :param str regalgpath: Path to registration algorithm executable.
        :param str regresultpath: Path where save registration results.
        :param list inputdims: List of input peak list dimension labels.
        :param list rootdims: List of root peak list dimension labels.
        :return: Results dict.
        :rtype: dict
        """
        # create command-line arguments for registration algorithm executable
        args = [plformat] + ["-noi"] + ["-save", regresultpath] + ["-dim"] + inputdims + [":"] + rootdims
        result = registration.run_registration(regalgpath, inputlistpath, rootlistpath, regresultpath, *args)
        return result

    def group(self, maxstep=10, maxregstep=2):
        """Group peaks from peak list into spin systems.

        :param int maxstep: Maximum number of iterations to perform.
        :param int maxregstep: Maximum number of iterations registration algorithm is allowed to run.
        :return:
        :rtype:
        """
        peaklistpath = self.peaklistpath
        peaklist = self.parse(self.peaklistpath, self.spectrumtype, self.dimlabels, self.plformat)
        peaklistfname = os.path.basename(self.peaklistpath)

        ids = (str(i) for i in range(1, 11))
        rdims = [label if label in self.rootdims else next(ids) for label in self.dimlabels]
        idims = [label if label in self.rootdims else next(ids) for label in self.dimlabels]

        initialreg_stds = {}
        registration_stds = {}
        previousrun_stds = {}
        currentrun_stds = {}
        usingreg = True
        istep = 0

        for i in range(0, maxstep):
            regresultpath = os.path.join(self.reg_result_dir, peaklistfname + str(istep) + ".json")
            groupresultpath = os.path.join(self.group_result_dir, peaklistfname + str(istep) + ".json")

            if istep >= maxregstep:
                usingreg = False

            if usingreg:
                regresult = self.calculate_registration(self.regalgpath, peaklistpath, peaklistpath,
                                                        self.formats[peaklist.plformat], regresultpath,
                                                        rdims, idims)

            if istep == 0:
                if not self._is_registered(regresult):
                    print("Peak list cannot be registered.")
                    break
                else:
                    initialreg_stds.update(regresult["FullSTD"])
                    previousrun_stds = dict(initialreg_stds)
                    registration_stds = {dimlabel: [] for dimlabel in initialreg_stds.keys()}

            if regresult and usingreg:
                peaklist_stds = regresult["FullSTD"]
                for dimlabel, std in peaklist_stds.items():
                    registration_stds[dimlabel].append(std)

                # do not let stds drop between iterations
                currentrun_stds = {dimlabel: max([previousrun_stds[dimlabel], peaklist_stds[dimlabel]])
                                   for dimlabel in previousrun_stds.keys()}

            else:  # (not regresult) or (regresult and not usingreg)
                usingreg = False
                currentrun_stds = {dimlabel: std+min(filter(None, registration_stds[dimlabel]))
                                   for dimlabel, std in currentrun_stds.items()}

            if not self._is_stds_within_range(currentrun_stds):
                break

            dbc = grouping.DBSCAN(datapath=self.peaklistpath, minpts=2)
            dbc.dbscan(data=peaklist, stds=currentrun_stds)
            dbc.print_clusters()

            with open(groupresultpath, "w") as outfile:
                dbc.write(outfile)

            previousrun_stds = dict(currentrun_stds)

            unclustered_peaks = dbc.get_noise()
            if unclustered_peaks:
                peaklist = pe.PeakList.fromlist(peaks=unclustered_peaks, spectrumtype=self.spectrumtype,
                                                dimlabels=self.dimlabels, plformat="json")
                peaklistpath = os.path.join(self.pl_dir, peaklistfname + str(istep) + ".json")
                with open(peaklistpath, "w") as outfile:
                    peaklist.write(filehandle=outfile, plformat="json")
            istep += 1

    def _is_registered(self, regresult):
        """Check if a peak list can be registered initially, i.e. is suitable for
        grouping into spin systems.

        :param dict regresult: Result of the first run of registration algorithm.
        :return: Peak list registered (True) or not registered (False).
        :rtype: bool
        """
        if not regresult:
            return False
        else:
            stds = regresult["FullSTD"]
            if not self._is_stds_within_range(stds):
                return False
        return True

    @staticmethod
    def _is_stds_within_range(stds, cutoff=0.5):
        """Check if each dimensions stds are within range, if the stds are too large they
        cannot be used for grouping.

        :param dict stds: Dictionary of stds for each dimension.
        :return: All dimension stds are valid (True) or not (False).
        :rtype: bool
        """
        if any([True for std in stds.values() if std >= cutoff]):
            return False
        return True

    @staticmethod
    def _load_json(filepath):
        """Load JSON file from a given file path.

        :param str filepath: Path to JSON file.
        :return: Dictionary represenation of JSON file.
        :rtype: dict
        """
        try:
            with open(filepath, "r") as infile:
                data = json.load(infile)
            return data
        except ValueError:
            raise ValueError('"{}" is not a valid JSON file.'.format(filepath))

if __name__ == "__main__":
    # python3 ssc.py ../datasets/jr19_hncocacb.pks

    script = sys.argv.pop(0)
    plpath = sys.argv.pop(0)
    stype = "HNcoCACB"
    labels = ["HN", "N", "CA/CB-1"]
    plfmt = "autoassign"
    rdims = ["HN", "N"]
    crspath = "./bin/calculate_registration"

    css = SpinSystemCreator(peaklistpath=plpath,
                            plformat=plfmt,
                            spectrumtype=stype,
                            dimlabels=labels,
                            rootdims=rdims,
                            regalgpath=crspath)
    css.group()
