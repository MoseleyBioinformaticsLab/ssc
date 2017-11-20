#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.peaklistparsers
~~~~~~~~~~~~~~~~~~~

This module has classes to parse NMR peak lists in different
    formats: Sparky, AutoAssign, JSON.
"""

import re
import json

from . import physicalentities as pe


def parse(peaklist_path, spectrum_type, labels, plformat):
    """Parse single peak list.

    :param str peaklist_path: Path to the peak list file.
    :param str spectrum_type: Type of the NMR experiment.
    :param list labels: List of dimension labels.
    :param str plformat: Peak list format.
    :return: Peak list.
    :rtype: :class:`~ssc.physicalentities.PeakList`
    """
    if plformat == "sparky":
        peaklist = SparkyPeakListParser.parse(peaklist_path, spectrum_type, labels, plformat)
    elif plformat == "autoassign":
        peaklist = AutoAssignPeakListParser.parse(peaklist_path, spectrum_type, labels, plformat)
    elif plformat == "json":
        peaklist = JSONPeakListParser.parse(peaklist_path, spectrum_type, labels, plformat)
    elif plformat == "cstable":
        peaklist = CSTablePeakListParser.parse(peaklist_path, spectrum_type, labels, plformat)
    else:
        raise TypeError('Unknown peak list format: "{}"'.format(plformat))

    filtered_peaklist = pe.PeakFilter.filterlist(peaklist=peaklist,
                                                 filters=[pe.ChemShiftPeakFilter({"CA": {"min": 35, "max": 75}}),
                                                          pe.ChemShiftPeakFilter({"N": {"min": 90, "max": 140}}),
                                                          pe.ChemShiftPeakFilter({"HN": {"min": 0, "max": 20}})])
    return filtered_peaklist


class PeakListParser(object):
    """Abstract PeakListParser class."""

    @classmethod
    def parse(cls, filepath, spectrum_type, labels, plformat):
        """Parse peak list content.

        :param str filepath: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of peak dimension labels.
        :param str plformat: Peak list format.
        """
        raise NotImplementedError("Subclass must implement abstract method.")

    @classmethod
    def write(cls, filehandle, peaklist, plformat):
        """Write peak list into file.

        :param filehandle: Peak list file pointer.
        :type filehandle: :py:class:`io.TextIOWrapper`
        :param peaklist: Peak list object.
        :type peaklist: :class:`~ssc.physicalentities.PeakList`
        :param str plformat: Peak list format.
        :return: None
        :rtype: None
        """
        try:
            if plformat == "json":
                json_str = cls.to_json(peaklist)
                filehandle.write(json_str)
            elif plformat == "autoassign":
                autoassign_str = cls.to_autoassign(peaklist)
                filehandle.write(autoassign_str)
            elif plformat == "sparky":
                sparky_str = cls.to_sparky(peaklist)
                filehandle.write(sparky_str)
            else:
                filehandle.close()
                raise TypeError("Unknown file format.")
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')

    @staticmethod
    def to_sparky(peaklist):
        """Convert peak list into a Sparky formatted string.

        :param peaklist: Peak list object.
        :type peaklist: :class:`~ssc.physicalentities.PeakList`
        :return: Sparky formatted string.
        :rtype: :py:class:`str`
        """
        sparky_str = "Assignment\t{}\n".format("".join(["w{}\t".format(i + 1) for i in range(0, len(peaklist.labels))]))
        for peak in peaklist:
            sparky_str += "-".join(peak.assignments_list) + "\t" + "\t".join(
                [str(cs) for cs in peak.chem_shifts_list]) + \
                          "\n"
        return sparky_str

    @staticmethod
    def to_autoassign(peaklist):
        """Convert peak list into a AutoAssign formatted string.

        :param peaklist: Peak list object.
        :type peaklist: :class:`~ssc.physicalentities.PeakList`
        :return: AutoAssign formatted string.
        :rtype: :py:class:`str`
        """
        autoassign_str = "Index\t{}Intensity\tWorkbook\n".format(
            "".join(["{}Dim\t".format(i + 1) for i in range(0, len(peaklist.labels))]))

        for index, peak in enumerate(peaklist):
            autoassign_str += str(index + 1) + "\t" + "\t".join([str(cs) for cs in peak.chem_shifts_list]) + "\t0\t" + \
                              peaklist.spectrum_type + "\n"
        return autoassign_str

    @staticmethod
    def to_json(peaklist):
        """Convert peak list into a JSON formatted string.

        :param peaklist: Peak list object.
        :type peaklist: :class:`~ssc.physicalentities.PeakList`
        :return: JSON formatted string.
        :rtype: :py:class:`str`
        """
        pl = [{"Assignment": peak.assignments_list,
               "Dimensions": peak.chem_shifts_list,
               "DataHeight": peak.extra_attr} for peak in peaklist]
        json_str = json.dumps(pl)
        return json_str


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
    def parse(cls, filepath, spectrum_type, labels, plformat):
        """Parse peak list in sparky format.

        :param str filepath: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of peak dimension labels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, 'r') as sparkyfile:
            peaklist = pe.PeakList(filepath, spectrum_type, labels, plformat)

            for line in sparkyfile:
                peak_attr = [float(x) if SparkyPeakListParser.numberPattern.match(x) else x for x in line.split()]
                if peak_attr == [] or not SparkyPeakListParser.assignmentPattern.match(peak_attr[0]):
                    continue

                assignment = re.split("[-]", peak_attr.pop(0))
                peaklist.append(pe.Peak(labels, assignment, peak_attr, peaklist))

        return peaklist


class AutoAssignPeakListParser(PeakListParser):
    """Concrete PeakListParser class that parses peak lists in autoassign format."""

    @classmethod
    def parse(cls, filepath, spectrum_type, labels, plformat):
        """Parse peak list in autoassign format.

        :param str filepath: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of peak dimension labels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, 'r') as autoassign_file:
            peaklist = pe.PeakList(filepath, spectrum_type, labels, plformat)

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
                peaklist.append(pe.Peak(labels, assignment, peak, peaklist))
        return peaklist


class JSONPeakListParser(PeakListParser):
    """Concrete PeakListParser class that parses peak lists in JSON format."""

    @classmethod
    def parse(cls, filepath, spectrum_type, labels, plformat):
        """Parse peak list in JSON format.

        :param str filepath: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of peak dimension labels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, "r") as json_file:
            data = json.load(json_file)
            peaklist = pe.PeakList(filepath, spectrum_type, labels, plformat)

            for peak in data:
                peaklist.append(pe.Peak(labels, peak["Assignment"], peak["Dimensions"], peaklist))
        return peaklist


class CSTablePeakListParser(PeakListParser):
    """Concrete PeakListParser class that parses peak list that consists of
    only columns of chemical shift values."""

    @classmethod
    def parse(cls, filepath, spectrum_type, labels, plformat):
        """Parse peak list in CSTable format.

        :param str filepath: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of peak dimension labels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        with open(filepath, "r") as cstable_file:
            peaklist = pe.PeakList(filepath, spectrum_type, labels, plformat)

            for line in cstable_file:
                if not line or line.startswith('#'):
                    continue

                peak = line.split()
                peak = [float(dim) for dim in peak]
                num_dims = len(peak)
                assignment = ["?" for _ in range(num_dims)]
                peaklist.append(pe.Peak(labels, assignment, peak, peaklist))

        return peaklist
