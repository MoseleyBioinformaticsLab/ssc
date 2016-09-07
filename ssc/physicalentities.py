#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
physicalentities.py

This module has classes to represent part of PhysicalEntities:
'PeakList', 'Peak', 'Dimension', 'Resonance'.

"""

import json
import pandas as pd


class PeakList(list):
    """Peak list container class."""

    formats = {"sparky": ".txt",
               "autoassign": ".pks",
               "json": ".json"}

    def __init__(self, filepath, spectrumtype, dimlabels, plformat):
        """Peak list initializer.

        :param str filepath: Path to the peak list file.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension dimlabels.
        :param str plformat: Peak list format.
        """
        super().__init__()
        self.filepath = filepath
        self.spectrumtype = spectrumtype
        self.dimlabels = dimlabels
        self.plformat = plformat

    @classmethod
    def fromlist(cls, peaks, spectrumtype, dimlabels, plformat):
        """Construct peak list object from a list of peaks.

        :param list peaks: List of peaks.
        :param str spectrumtype: Type of the NMR experiment.
        :param list dimlabels: List of dimension dimlabels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        peaklist = cls(filepath=None, spectrumtype=spectrumtype, dimlabels=dimlabels, plformat=plformat)
        for peak in peaks:
            peaklist.append(peak)
        return peaklist

    def write(self, filehandle, plformat):
        """Save :class:`~ssc.physicalentities.PeakList` into file.

        :param filehandle: Peak list file pointer.
        :param str plformat: Peak list format.
        :return: None
        :rtype: None
        """
        try:
            if plformat == "json":
                json_str = self._to_json()
                filehandle.write(json_str)
            elif plformat == "autoassign":
                autoassign_str = self._to_autoassign()
                filehandle.write(autoassign_str)
            elif plformat == "sparky":
                sparky_str = self._to_sparky()
                filehandle.write(sparky_str)
            else:
                filehandle.close()
                raise TypeError("Unknown file format.")
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')

    def _to_sparky(self):
        """Save :class:`~ssc.physicalentities.PeakList` into Sparky formatted string.

        :return: Sparky formatted string.
        :rtype: str
        """
        sparky_str = "Assignment\t{}\n".format("".join(["w{}\t".format(i+1) for i in range(0, len(self.dimlabels))]))
        for peak in self:
            sparky_str += "-".join(peak.assignmentslist) + "\t" + "\t".join([str(cs) for cs in peak.chemshiftslist]) + \
                          "\n"
        return sparky_str

    def _to_autoassign(self):
        """Save :class:`~ssc.physicalentities.PeakList` into AutoAssign formatted string.

        :return: AutoAssign formatted string.
        :rtype: str
        """
        autoassign_str = "Index\t{}Intensity\tWorkbook\n".format(
            "".join(["{}Dim\t".format(i+1) for i in range(0, len(self.dimlabels))]))

        for index, peak in enumerate(self):
            autoassign_str += str(index+1) + "\t" + "\t".join([str(cs) for cs in peak.chemshiftslist]) + "\t0\t" + \
                              self.spectrumtype + "\n"
        return autoassign_str

    def _to_json(self):
        """Save :class:`~ssc.physicalentities.PeakList` into JSON formatted string.

        :return: JSON formatted string.
        :rtype: str
        """
        peaklist = [{"Assignment": peak.assignmentslist,
                     "Dimensions": peak.chemshiftslist,
                     "DataHeight": peak.extra_attr} for peak in self]
        json_str = json.dumps(peaklist, sort_keys=False, indent=4)
        return json_str

    @property
    def peaklistdf(self):
        """DataFrame representation of a peak list.

        :return: DataFrame representation of a peak list.
        :rtype: :class:`~pandas.DataFrame`
        """
        return pd.DataFrame([peak.chemshiftslist for peak in self], columns=self.dimlabels)


class Peak(list):
    """Peak container class."""
    
    def __init__(self, labels, assignment, peak_attr, owner):
        """Peak initializer.

        :param list labels: List of dimension labels.
        :param list assignment: List of dimension assignments.
        :param list peak_attr: List of peak attributes.
        :param :class:`~ssc.physicalentities.PeakList` owner: Peak list object where peak belongs.
        """
        super().__init__()
        self.labels = labels
        self.assignment = assignment
        self.peak_attr = peak_attr
        self.owner = owner

        for idx, label in enumerate(labels):
            super().append(Dimension(idx + 1, label, assignment[idx], peak_attr[idx]))

        self.extra_attr = peak_attr[len(assignment):]

    @property
    def chemshiftslist(self):
        """List of chemical shifts.

        :return: List of chemical shifts.
        :rtype: list
        """
        return [dim.chemshift for dim in self]

    @property
    def assignmentslist(self):
        """List of assignments.

        :return: List of assignments.
        :rtype: list
        """
        return [dim.assignment for dim in self]

    @property
    def chemshiftsdict(self):
        """Dictionary of chemical shifts.

        :return: Dictionary of label-chemical shift key-value pairs.
        :rtype: dict
        """
        return {label: chemshift for label, chemshift in zip(self.labels, self.chemshiftslist)}

    @property
    def assignmentsdict(self):
        """Dictionary of assignments.

        :return: Dictionary of label-assignment key-value pairs.
        :rtype: dict
        """
        return {label: assignment for label, assignment in zip(self.labels, self.assignmentslist)}


class Dimension(object):
    """Class that represents each dimension of a peak."""

    def __init__(self, dimid, label, assignment, chemshift):
        """Dimension initializer.

        :param int dimid: Dimension index.
        :param str label: Dimension label.
        :param str assignment: Dimension assignment.
        :param float chemshift: Dimension chemical shift value.
        """
        self.dimid = dimid
        self.label = label
        self.assignment = assignment
        self.chemshift = chemshift

    def is_assigned(self):
        """Test if dimension is assigned.

        :return: Assigned (True), not assigned (False).
        :rtype: bool
        """
        return self.assignment != '?' and self.assignment != ''

    def __str__(self):
        """String representation of dimension.

        :return: String representation of dimension.
        :rtype: str
        """
        return "{dimid}{assignment}{chemshift}".format(**self.__dict__)

    def __repr__(self):
        """String representation of dimension.

        :return: String representation of dimension.
        :rtype: str
        """
        return str(self.chemshift)


class Resonance(object):
    """Resonance class"""
    def __init__(self, chemshift):
        self.chemshift = chemshift


class PeakListFilter(object):

    @staticmethod
    def filter(peaklist, filters):
        """Apply multiple filters to a peak list.

        :param peaklist: Peak list to be filtered.
        :param dict filters: List of peak filters.
        :type peaklist: :class:`~ssc.physicalentities.PeakList`
        :return: Filtered peak list.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        peaks = list(filter(lambda peak: all(f(peak, f_param) for f, f_param in filters.items()), peaklist))
        return PeakList.fromlist(peaks, peaklist.spectrumtype, peaklist.dimlabels, peaklist.plformat)

    @staticmethod
    def peak_chemshift_filter(peak, dim_range):
        """Filter peak(s) based on min and max chemical shift value for a given dimension.

        :param peak: Single peak within peak list.

        :param dict dim_range: Dictionary specifying min and max values for dim label.
        :return: If peak passes chemical shift filter (True) or not (False).
        :rtype: bool
        """
        return all([dim_range[label]["min"] <= peak.chemshiftsdict[label] <= dim_range[label]["max"]
                    for label in dim_range.keys() if label in peak.labels])
