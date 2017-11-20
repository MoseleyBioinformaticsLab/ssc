#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.grouping
~~~~~~~~~~~~

This module implements algorithm for grouping peaks into proto-spin systems
using modified DBSCAN clustering algorithm: it uses registration algorithm
to calculate stds necessary to group peaks that belong to the same spin system
into clusters and uses chi-square probability formalism to decide if peak
belongs to the spin system cluster.
"""

import json
from math import sqrt

from scipy.stats import chi2


class Cluster(object):
    """A cluster container for spin systems, i.e. groups of peaks (points)
    that are clustered together and belong to the same spin system."""

    def __init__(self, label, stds):
        """Initialize cluster.

        :param int label: Cluster label (id).
        :param dict stds: Stds that were used to group peaks into cluster.
        """
        self.label = label
        self.stds = stds
        self.members = []

    def add_member(self, member):
        """Add member into cluster.

        :param member: Represents single peak within PeakList object.
        :type: :class:`~ssc.physicalentities.Peak`
        :return: None.
        :rtype: :py:obj:`None`
        """
        self.members.append(member)

    def coordinates(self, label):
        """Construct list of coordinates for cluster for a specified dimension label.

        :param str label: Peak dimension label.
        :return: List of coordinates.
        :rtype: :py:class:`list`
        """
        return [peak.chem_shifts_dict[label] for peak in self.members]

    def assignments(self):
        """Construct list of assignments for each peak member in a cluster.

        :return: List of assignments.
        :rtype: :py:class:`list`
        """
        return [peak.assignments_list for peak in self.members]

    def __contains__(self, member):
        """Test if member is inside members list.

        :param member: Represents single peak within PeakList object.
        :type member: :class:`~ssc.physicalentities.Peak`
        :return: Peak is inside members list (True) or or not (False).
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        return member in self.members


class DBSCAN(object):
    """Modified DBSCAN clustering algorithm that uses experimental NMR peak list
    and standard deviations for comparable peak list dimensions in order to
    identify peaks that belong to the same spin system (clustered together).
    """
    global_clusters = {}

    def __init__(self, data_path, min_pts):
        """Initialize clustering algorithm.

        :param str data_path: Path to peak list.
        :param int min_pts: Minimum number of points that can form a cluster.
        """
        self.global_clusters.setdefault(data_path, {})
        self.count = max(self.global_clusters[data_path].keys(), default=0)
        self.visited = []
        self.clusters = {}
        self.min_pts = min_pts
        self.data_path = data_path
        self.noise = Cluster(label=-1, stds={})

    def dbscan(self, data, stds):
        """Run modified DBSCAN algorithm.

        :param data: Represents experimental NMR peak list.
        :type data: :class:`~ssc.physicalentities.PeakList`
        :param dict stds: Stds used to group :class:`~ssc.physicalentities.Peak` into cluster.
        :return: List of clusters.
        :rtype: :py:class:`list`
        """
        for point in data:
            if point not in self.visited:
                self.visited.append(point)
                region_pts = self.region_query(point, data, stds)
                if len(region_pts) < self.min_pts:
                    self.noise.add_member(point)
                else:
                    self.count += 1
                    cluster = Cluster(label=self.count, stds=stds)
                    self.expand_cluster(point, region_pts, cluster, data, stds)

        self.clusters[-1] = self.noise
        self.global_clusters[self.data_path][-1] = self.noise
        return self.clusters

    @staticmethod
    def region_query(point, data, stds, prob=0.0001):
        """Find a region around point that satisfies given probability criteria.
        We are using chi-square probability to decide if a point belongs to a
        region or not.

        :param point: Represents single peak within PeakList object.
        :param data: Represents PeakList object.
        :param dict stds: Stds used to group points into cluster.
        :param float prob: Probability.
        :type point: :class:`~ssc.physicalentities.Peak`
        :type data: :class:`~ssc.physicalentities.PeakList`
        :return: List of neighbor points.
        :rtype: :py:class:`list`
        """
        region = []
        df = len(stds)  # degrees of freedom
        chi2dist_cutoff = sqrt(chi2.isf(prob, df))

        # prevent 0-division
        stds = {label: stds[label] if stds[label] > 0 else 0.0000001 for label in stds}

        for otherpoint in data:
            dist = 0
            for idx, label, in enumerate(data.labels):
                if label in stds:
                    dist += ((otherpoint[idx].chem_shift - point[idx].chem_shift) / stds[label]) ** 2
            normalized_dist = sqrt(dist)
            if normalized_dist <= chi2dist_cutoff:
                region.append(otherpoint)
        return region

    def expand_cluster(self, point, neighbor_pts, cluster, data, stds):
        """Expand cluster by examining neighborhood region and add_member it to list of clusters.

        :param point: Represents single peak within peak list.
        :param list neighbor_pts: List of neighbor points.
        :param cluster: Represents particular cluster of points.
        :param data: Represents PeakList object.
        :param dict stds: Stds used to group points into cluster.
        :type point: :class:`~ssc.physicalentities.Peak`
        :type cluster: :class:`~ssc.grouping.Cluster`
        :return: None
        :rtype: :py:obj:`None`
        """
        cluster.add_member(point)
        for n_pt in neighbor_pts:
            if n_pt not in self.visited:
                self.visited.append(n_pt)
                new_neighbor_pts = self.region_query(n_pt, data, stds)
                if len(new_neighbor_pts) >= self.min_pts:
                    for nn_pt in new_neighbor_pts:
                        if nn_pt not in neighbor_pts:
                            neighbor_pts.append(nn_pt)

            for clstr in self.clusters.values():
                if n_pt not in clstr:
                    if n_pt not in cluster:
                        cluster.add_member(n_pt)

            # if n_pt not yet member of any cluster, add n_pt to cluster
            if len(self.clusters) == 0:
                if n_pt not in cluster:
                    cluster.add_member(n_pt)

        self.clusters[self.count] = cluster
        self.global_clusters[self.data_path][self.count] = cluster

    def write(self, filehandle, outformat="json"):
        """Write clustering results into file.

        :param filehandle: Results file pointer.
        :param str outformat: Rusults file format to use.
        :return: None.
        :rtype: :py:obj:`None`
        """
        try:
            if outformat == "json":
                json_str = self.to_json()
                filehandle.write(json_str)
            else:
                filehandle.close()
                raise TypeError("Unknown file format.")
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')

    def to_json(self):
        """Save list of clusters into JSON formatted string.

        :return: List of clusters as JSON formatted string.
        :rtype: :py:class:`str`
        """
        clusters = []
        for cluster in self.global_clusters[self.data_path].values():
            cluster_dict = {"label": cluster.label,
                            "peaks": [],
                            "stds": cluster.stds}
            for idx, peak in enumerate(cluster.members):
                peak_dict = {"dimensions": peak.chem_shifts_list,
                             "assignment": peak.assignments_list,
                             "index": idx}
                cluster_dict["peaks"].append(peak_dict)
            clusters.append(cluster_dict)

        json_str = json.dumps(clusters, sort_keys=False, indent=4)
        return json_str
