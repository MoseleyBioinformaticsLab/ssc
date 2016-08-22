#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
physicalentities.py

This module has classes to represent part of PhysicalEntities group, e.g.
'PeakList', 'Peak', 'Dimension', 'Resonance'.

"""

import abc
import json


class PeakList(list):
    """Class representing experimental peak list."""

    formats = {"sparky": ".txt",
               "autoassign": ".pks",
               "json": ".json"}

    def __init__(self, filepath, spectrumtype, dimlabels, plformat):
        """Initialization method.

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

        # self.notUsedPeaks = []

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
            peaklist.add_peak(peak)
        return peaklist

    def add_peak(self, peak):
        """Add peak object into peak list object.

        :param peak: Peak object.
        :type peak: :class:`~ssc.physicalentities.Peak`
        :return: None
        :rtype: None
        """
        peak.owner(self)
        super().append(peak)

    def write(self, filehandle, plformat):
        """Save :class:`~ssc.physicalentities.PeakList` data into file object.

        :param filehandle:
        :param plformat:
        :return:
        """
        try:
            if plformat is "json":
                json_str = self._to_json()
                filehandle.write(json_str)
            elif plformat is "autoassign":
                autoassign_str = self._to_autoassign()
                filehandle.write(autoassign_str)
            elif plformat is "sparky":
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
                     "DataHeight": peak.extraattr} for peak in self]
        json_str = json.dumps(peaklist, sort_keys=False, indent=4)
        return json_str


class Peak(list):
    """Peak class."""
    
    def __init__(self, assignment, labels, peakattr):
        """Peak initializer.

        :param assignment:
        :param labels:
        :param peakattr:
        """
        super().__init__()
        self.assignment = assignment
        self.labels = labels
        self.peakAttr = peakattr

        for dimension, label in zip(range(0, len(assignment)), labels):
            super().append(Dimension(dimension + 1, label, assignment[dimension], peakattr[dimension]))

        self.extraattr = peakattr[len(assignment):]

    def owner(self, peakList=None):
        """Method that sets the 'PeakList' object as owner of 'Peak' object."""
        if peakList != None:
            self.peakList = peakList
        return self.peakList

    @property
    def chemshiftslist(self):
        return [dim.chemshift for dim in self]

    @property
    def assignmentslist(self):
        return [dim.assignment for dim in self]


class Dimension(object):
    """Dimension class is used by Peak object to represent dimensions within a 
    given peak."""

    def __init__(self, dimid, label, assignment, chemshift):
        """Dimension initializer.

        :param dimid:
        :param label:
        :param assignment:
        :param chemshift:
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
        """String representation of :class:`~ssc.physicalentities.Dimension`."""
        return "{dimID}{dimAssignment}{chemShift}".format(**self.__dict__)

    def __repr__(self):
        return str(self.chemshift)


class Resonance(object):
    """Resonance class"""
    def __init__(self, chemshift):
        self.chemshift = chemshift


class PeakListFilter(metaclass=abc.ABCMeta):

    def __init__(self, peakfilter):
        self.peakfilter = peakfilter

    @abc.abstractmethod
    def filteredPeakList(self, peakList):
        """Filter peaks from peak list

        :param peakList: Peak list
        :type peakList: sass.pe.PeakList
        :return: Filetered peak list
        :rtype: sass.pe.PeakList
        """
        raise NotImplementedError("Subclass must implement abstract method.")

class DimRangePeakListFilter(PeakListFilter):

    def filteredPeakList(self, peakList):
        """Filter peaks from peak list

        :param peakList: Peak list
        :type peakList: sass.pe.PeakList
        :return: Filetered peak list
        :rtype: sass.pe.PeakList
        """
        notUsedPeaks = []

        for peak in peakList:
            withinRange = self.checkDimRange(peak)
            if withinRange:
                continue
            elif not withinRange:
                notUsedPeaks.append(peak)

        for peak in notUsedPeaks:
            if peak in peakList:
                peakList.remove(peak)

        return peakList

    def checkDimRange(self, peak):
        """Check if peak's comparable dimensions are within range

        :param peak: Peak object
        :type peak: sass.pe.Peak
        :param dict filter: Dimension range filter
        :return: If peak dimensions are within allowed dim range
        :rtype: bool
        """
        result = []
        peak = [dim for dim in peak if dim.dimLabel in self.peakfilter]
        for dim in peak:
            if (self.peakfilter[dim.dimLabel]['min'] <= dim.chemShift <= self.peakfilter[dim.dimLabel]['max']):
                result.append(True)
            else:
                result.append(False)
        return all(result)
