#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.peaklistparsers
~~~~~~~~~~~~~~~~~~~

This module has classes to parse NMR peak lists in different formats: Sparky, AutoAssign, JSON.
"""

import abc
import re
import json

import physicalentities as pe


class PeakListParser(metaclass=abc.ABCMeta):
    """Abstract PeakListParser class."""

    @abc.abstractmethod
    def parse(self, filepath, spectrumtype, dimlabels, plformat):
        """Parse peak list content.

        :param str filepath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension dimlabels.
        :param str plformat: Peak list format.
        """
        raise NotImplementedError("Subclass must implement abstract method.")


class SparkyPeakListParser(PeakListParser):
    """Concrete PeakListParser class that parses peak lists in sparky format."""

    numberPattern = re.compile(r"""
                                \d*             # integral part
                                [.]?            # decimal point
                                \d+             # fractional part
                                """, re.VERBOSE)

    assignmentPattern = re.compile(r"""
                                    [\w\?]+         # 1st dimension assignment
                                    ([-][\w\?+])+   # 2nd, 3rd or 4th dimension, separated by "-"
                                    """, re.VERBOSE)

    @classmethod
    def parse(cls, filepath, spectrumtype, dimlabels, plformat):
        """Parse peak list in sparky format.

        :param str filepath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension dimlabels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, 'r') as sparkyfile:
            peaklist = pe.PeakList(filepath, spectrumtype, dimlabels, plformat)

            for line in sparkyfile:
                peak_attr = [float(x) if SparkyPeakListParser.numberPattern.match(x) else x for x in line.split()]
                if peak_attr == [] or not SparkyPeakListParser.assignmentPattern.match(peak_attr[0]):
                    continue

                assignment = re.split("[-]", peak_attr.pop(0))
                peaklist.append(pe.Peak(dimlabels, assignment, peak_attr, peaklist))
        
        return peaklist


class AutoAssignPeakListParser(PeakListParser):
    """Concrete PeakListParser class that parses peak lists in autoassign format."""

    @classmethod
    def parse(cls, filepath, spectrumtype, dimlabels, plformat):
        """Parse peak list in autoassign format.

        :param str filepath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension dimlabels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, 'r') as autoassign_file:
            peaklist = pe.PeakList(filepath, spectrumtype, dimlabels, plformat)

            for line in autoassign_file:
                if line.startswith('#'):
                    continue
                if line.startswith('*'):
                    break

                peak = line.split()
                index = peak.pop(0)
                workbook = peak.pop()
                intensity = peak.pop()
                num_dims = len(peak)
                assignment = ['?' for i in range(num_dims)]
                peak = [float(dim) for dim in peak]
                peaklist.append(pe.Peak(dimlabels, assignment, peak, peaklist))

        return peaklist


class JSONPeakListParser(PeakListParser):
    """Concrete PeakListParser class that parses peak lists in JSON format."""

    @classmethod
    def parse(cls, filepath, spectrumtype, dimlabels, plformat):
        """Parse peak list in JSON format.

        :param str filepath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension dimlabels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, "r") as json_file:
            data = json.load(json_file)
            peaklist = pe.PeakList(filepath, spectrumtype, dimlabels, plformat)

            for peak in data:
                peaklist.append(pe.Peak(dimlabels, peak["Assignment"], peak["Dimensions"], peaklist))

        return peaklist
