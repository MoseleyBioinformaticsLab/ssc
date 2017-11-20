#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.physicalentities
~~~~~~~~~~~~~~~~~~~~

This module has classes to represent part of physical entities:
peak list, peak, dimension, resonance.
"""

import pandas as pd


class PeakList(list):
    """Peak list container class."""

    formats = {"sparky": ".txt",
               "autoassign": ".pks",
               "json": ".json"}

    def __init__(self, filepath, spectrum_type, labels, plformat):
        """Peak list initializer.

        :param str filepath: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of dimension labels.
        :param str plformat: Peak list format.
        """
        super().__init__()
        self.filepath = filepath
        self.spectrum_type = spectrum_type
        self.labels = labels
        self.plformat = plformat

    @classmethod
    def fromlist(cls, peaks, spectrum_type, labels, plformat):
        """Construct peak list object from a list of peaks.

        :param list peaks: List of peaks.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of dimension labels.
        :param str plformat: Peak list format.
        :return: Peak list object.
        :rtype: :class:`~ssc.physicalentities.PeakList`
        """
        peaklist = cls(filepath=None, spectrum_type=spectrum_type, labels=labels, plformat=plformat)
        for peak in peaks:
            peaklist.append(peak)
        return peaklist

    @property
    def peaklistdf(self):
        """DataFrame representation of a peak list.

        :return: DataFrame representation of a peak list.
        :rtype: :class:`~pandas.DataFrame`
        """
        return pd.DataFrame([peak.chem_shifts_list for peak in self], columns=self.labels)


class Peak(list):
    """Peak container class."""

    def __init__(self, labels, assignment, peak_attr, owner):
        """Peak initializer.

        :param list labels: List of dimension labels.
        :param list assignment: List of dimension assignments.
        :param list peak_attr: List of peak attributes.
        :param owner: Peak list object where peak belongs.
        :type owner: :class:`~ssc.physicalentities.PeakList`
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
    def chem_shifts_list(self):
        """List of chemical shifts.

        :return: List of chemical shifts.
        :rtype: :py:class:`list`
        """
        return [dim.chem_shift for dim in self]

    @property
    def assignments_list(self):
        """List of assignments.

        :return: List of assignments.
        :rtype: :py:class:`list`
        """
        return [dim.assignment for dim in self]

    @property
    def chem_shifts_dict(self):
        """Dictionary of chemical shifts.

        :return: Dictionary of label-chemical shift key-value pairs.
        :rtype: :py:class:`dict`
        """
        return {label: chemshift for label, chemshift in zip(self.labels, self.chem_shifts_list)}

    @property
    def assignments_dict(self):
        """Dictionary of assignments.

        :return: Dictionary of label-assignment key-value pairs.
        :rtype: :py:class:`dict`
        """
        return {label: assignment for label, assignment in zip(self.labels, self.assignments_list)}


class Dimension(object):
    """Class that represents each dimension of a peak."""

    def __init__(self, dim_id, label, assignment, chem_shift):
        """Dimension initializer.

        :param int dim_id: Dimension index.
        :param str label: Dimension label.
        :param str assignment: Dimension assignment.
        :param float chem_shift: Dimension chemical shift value.
        """
        self.dim_id = dim_id
        self.label = label
        self.assignment = assignment
        self.chem_shift = chem_shift

    def is_assigned(self):
        """Test if dimension is assigned.

        :return: Assigned (True), not assigned (False).
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        return self.assignment != "?" and self.assignment != ""

    def __str__(self):
        """String representation of dimension.

        :return: String representation of dimension.
        :rtype: :py:class:`str`
        """
        return "{dim_id}{assignment}{chem_shift}".format(**self.__dict__)

    def __repr__(self):
        """String representation of dimension.

        :return: String representation of dimension.
        :rtype: :py:class:`str`
        """
        return str(self.chem_shift)


class Resonance(object):
    """Resonance class"""

    def __init__(self, chem_shift):
        self.chem_shift = chem_shift


class PeakFilter(object):
    def __init__(self, filter_parameters):
        self.filter_parameters = filter_parameters

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
        return PeakList.fromlist(peaks, peaklist.spectrum_type, peaklist.labels, peaklist.plformat)


class ChemShiftPeakFilter(PeakFilter):
    def filter(self, peak):
        """Filter peak(s) based on min and max chemical shift value for a given dimension.

        :param peak: Single peak within peak list.
        :return: If peak passes chemical shift filter (True) or not (False).
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        return all([self.filter_parameters[label]["min"] <= peak.chem_shifts_dict[label] <=
                    self.filter_parameters[label]["max"]
                    for label in self.filter_parameters.keys() if label in peak.labels])
