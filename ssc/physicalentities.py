#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
physicalentities.py

This module has classes to represent part of PhysicalEntities:
'PeakList', 'Peak', 'Dimension', 'Resonance'.

"""

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


class PeakFilter(object):

    def __init__(self, filterparam):
        self.filterparam = filterparam

    @staticmethod
    def filterlist(peaklist, filters):
        """Apply multiple filters to a peak list.

        :param peaklist: Peak list to be filtered.
        :param list filters: List of peak list filters.
        :type peaklist: :class:`~ssc.physicalentities.PeakList`
        :return: Filtered peak list.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        peaks = list(filter(lambda peak: all(f.filter(peak) for f in filters), peaklist))
        return PeakList.fromlist(peaks, peaklist.spectrumtype, peaklist.dimlabels, peaklist.plformat)

class ChemShiftPeakFilter(PeakFilter):

    def filter(self, peak):
        """Filter peak(s) based on min and max chemical shift value for a given dimension.

        :param peak: Single peak within peak list.
        :return: If peak passes chemical shift filter (True) or not (False).
        :rtype: bool
        """
        return all([self.filterparam[label]["min"] <= peak.chemshiftsdict[label] <= self.filterparam[label]["max"]
                    for label in self.filterparam.keys() if label in peak.labels])
