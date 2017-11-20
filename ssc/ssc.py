#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.ssc
~~~~~~~

This module provides :class:`~ssc.ssc.SpinSystemCreator` class that groups peaks from
a single peak list into spin systems.
"""

import os
import json

from . import peaklistparsers as plp
from . import physicalentities as pe
from . import registration
from . import grouping


class SpinSystemCreator(object):
    """Perform single peak list registration and single peak list  grouping,
    i.e. group peaks that belong to the same spin system into clusters."""

    def __init__(self, peaklist_path, plformat, spectrum_type, labels, root_dims, regalg_path, grouping_result_path=""):
        """Initialize spin system creator.

        :param str peaklist_path: Path to the peak list file.
        :param str spectrum_type: Type of the NMR experiment.
        :param list labels: List of dimension labels for a given peak list.
        :param list root_dims: List of root dimension labels for a given peak list.
        :param str plformat: Peak list format.
        :param str regalg_path: Path to registration algorithm executable.
        :param str grouping_result_path: Path to where save grouping algorithm results.
        """
        self.peaklist_path = os.path.normpath(peaklist_path)
        self.spectrum_type = spectrum_type
        self.labels = labels
        self.root_dims = root_dims
        self.plformat = plformat
        self.regalg_path = os.path.normpath(regalg_path)

        peaklist_fname = os.path.basename(peaklist_path)
        grouping_result_fname = peaklist_fname + "_grouping_result.json"

        if not grouping_result_path:
            self.grouping_result_path = os.path.join(os.path.dirname(self.peaklist_path), grouping_result_fname)
        else:
            if not os.path.exists(grouping_result_path):
                os.mkdir(grouping_result_path)

            self.grouping_result_path = os.path.join(grouping_result_path, grouping_result_fname)

    def group(self, max_step=10, max_reg_step=2, view=False):
        """Group peaks from peak list into spin systems.

        :param int max_step: Maximum number of iterations grouping algorithm is allowed to run.
        :param int max_reg_step: Maximum number of iterations registration algorithm is allowed to run.
        :return: JSON string with final grouping results.
        :rtype: :py:class:`str`
        """
        peaklist = plp.parse(self.peaklist_path, self.spectrum_type, self.labels, self.plformat)
        peaklist_json_str = json.dumps(plp.PeakListParser.to_json(peaklist))

        ids = (str(i) for i in range(100, 200))
        rdims = [label if label in self.root_dims else next(ids) for label in self.labels]
        idims = [label if label in self.root_dims else next(ids) for label in self.labels]

        initial_registration_stds = {}
        registration_stds = {}
        previous_run_stds = {}
        current_run_stds = {}
        using_registration = True
        step = 0

        for i in range(0, max_step):
            if step >= max_reg_step:
                using_registration = False

            if using_registration:
                regresult = self.calculate_registration(self.regalg_path, peaklist_json_str, peaklist_json_str, idims, rdims)

            if step == 0:
                if not self._is_registered(regresult):
                    print("Peak list cannot be registered.")
                    break
                else:
                    initial_registration_stds.update(regresult["FullSTD"])
                    previous_run_stds = dict(initial_registration_stds)
                    registration_stds = {label: [] for label in initial_registration_stds.keys()}

            if regresult and using_registration:
                peaklist_stds = regresult["FullSTD"]
                for dimlabel, std in peaklist_stds.items():
                    registration_stds[dimlabel].append(std)

                # do not let stds drop between iterations
                current_run_stds = {label: max([previous_run_stds[label], peaklist_stds[label]]) for label in previous_run_stds.keys()}

            else:
                using_registration = False
                current_run_stds = {label: std+min(filter(None, registration_stds[label])) for label, std in current_run_stds.items()}

            if not self._is_stds_within_range(current_run_stds):
                break

            dbc = grouping.DBSCAN(data_path=self.peaklist_path, min_pts=2)
            dbc.dbscan(data=peaklist, stds=current_run_stds)
            previous_run_stds = dict(current_run_stds)
            unclustered_peaks = dbc.noise.members

            if unclustered_peaks:  # try to group until we still have unclustered peaks
                peaklist = pe.PeakList.fromlist(peaks=unclustered_peaks, spectrum_type=self.spectrum_type, labels=self.labels, plformat="json")
                peaklist_json_str = json.dumps(plp.PeakListParser.to_json(peaklist))
                step += 1
            else:
                break  # no peaks left to group

        with open(self.grouping_result_path, "w") as outfile:
            grouping_result_str = dbc.to_json()
            outfile.write(grouping_result_str)

        if view:
            print(grouping_result_str)

        return grouping_result_str

    @staticmethod
    def calculate_registration(regalg_path, input_peaklist, root_peaklist, input_dims, root_dims):
        """Execute registration algorithm binary.

        :param str input_peaklist: Path to input peak list.
        :param str root_peaklist: Path to root peak list.
        :param str regalg_path: Path to registration algorithm executable.
        :param list input_dims: List of input peak list dimension labels.
        :param list root_dims: List of root peak list dimension labels.
        :return: Dictionary containing registration results.
        :rtype: :py:class:`dict`
        """
        # create command-line arguments for registration algorithm executable
        args = ["--noi"] + ["--dim"] + input_dims + [":"] + root_dims
        result = registration.run_registration(regalg_path, input_peaklist, root_peaklist, *args)
        return result

    def _is_registered(self, registration_result):
        """Check if a peak list can be registered initially, i.e. is suitable for
        grouping into spin systems.

        :param dict registration_result: Result of the first run of the registration algorithm.
        :return: Peak list registered (True) or not registered (False).
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        if not registration_result:
            return False
        else:
            stds = registration_result["FullSTD"]
            if not self._is_stds_within_range(stds):
                return False
        return True

    @staticmethod
    def _is_stds_within_range(stds, cutoff=0.5):
        """Check if each dimensions stds are within range, if the stds are too large they
        cannot be used for grouping.

        :param dict stds: Dictionary of stds for each dimension.
        :return: All dimension stds are valid (True) or not (False).
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        if any([True for std in stds.values() if std >= cutoff]):
            return False
        return True
