#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.grouping
~~~~~~~~~~~~~~

This module implements algorithm for grouping peaks into proto-spin systems
using modified DBSCAN clustering algorithm.

"""

import json
from math import sqrt

from scipy.stats import chi2
import numpy as np


class Cluster(object):
    """A cluster container for spin systems, i.e. groups of peaks (points)
    that are clustered together and belong to the same spin system."""

    def __init__(self, label, stds):
        """Initialize Cluster.

        :param int label: Cluster label (id).
        :param dict stds: Stds that were used to group peaks into Cluster.
        """
        self.label = label
        self.stds = stds
        self.members = []
        self.x = []
        self.y = []
        self.desc = []
        self.overlapped = False

    def add_member(self, member):
        """Add member into Cluster.

        :param member: Represents single peak within PeakList object.
        :type: sass.pe.PhysicalEntities.Peak
        :return: None
        """
        self.members.append(member)
        self.x.append(member[0].chemshift)
        self.y.append(member[1].chemshift)
        self.desc.append(member.assignment)

    def __contains__(self, member):
        """Test if member is inside members_list.

        :param member: Represents single peak within PeakList object.
        :type member: sass.pe.PhysicalEntities.Peak
        :return: True or False
        :rtype: bool
        """
        return member in self.members

    # @property
    # def representative(self):
    #     """Representative point of a cluster by calculating mean across dimensions.
    #
    #     :return: Representative point.
    #     :rtype:
    #     """
    #     # points = np.array([[dim.chemshift for dim in point if dim.label in self.stds] for point in self.members])
    #     points = np.array([point.chemshiftslist for point in self.members])
    #     representative = np.mean(points, axis=0)
    #     return representative


class DBSCAN(object):
    """Modified DBSCAN clustering algorithm that uses experimental NMR peak list
    and standard deviations for comparable peak list dimensions in order to
    identify peaks that belong to the same spin system (clustered together).
    """

    global_clusters = {}

    def __init__(self, datapath, minpts):
        """Initialize clustering algorithm.

        :param str datapath: Path to initial peak list.
        :param int minpts: Minimum number of points that can form a cluster.
        """
        self.global_clusters.setdefault(datapath, {})
        self.count = max(self.global_clusters[datapath].keys(), default=0)
        self.visited = []
        self.clusters = {}
        self.minpts = minpts
        self.datapath = datapath
        self.noise = Cluster(label=-1, stds={})

    def dbscan(self, data, stds):
        """Run modified DBSCAN algorithm.

        :param data: Represents experimental NMR peak list.
        :type data: :class:`~ssc.physicalentities.PeakList`
        :param dict stds: Stds used to group :class:`~ssc.physicalentities.Peak` into cluster.
        :return: List of clusters.
        :rtype: list
        """
        for point in data:
            if point not in self.visited:
                self.visited.append(point)
                region_pts = self.region_query(point, data, stds)
                if len(region_pts) < self.minpts:
                    self.noise.add_member(point)
                else:
                    cluster = Cluster(label=str(self.count), stds=stds)
                    self.count += 1
                    self.expand_cluster(point, region_pts, cluster, data, stds)

        self.clusters[-1] = self.noise
        self.global_clusters[self.datapath][-1] = self.noise
        return self.clusters

    def region_query(self, point, data, stds, prob=0.0001):
        """Find a region around point that satisfies given probability criteria.
        We are using chi square probability to decide if a point belongs to a
        region or not.

        :param point: Represents single peak within PeakList object.
        :type point: :class:`~ssc.physicalentities.Peak`
        :param data: Represents PeakList object.
        :type data: :class:`~ssc.physicalentities.PeakList`
        :param dict stds: Stds used to group :class:`~ssc.physicalentities.Peak` into cluster.
        :param float prob: Probability.
        :return: List of neighbor points.
        :rtype: list
        """
        region = []
        df = len(stds)  # degrees of freedom
        chi2dist = sqrt(chi2.isf(prob, df))

        for otherpoint in data:
            dist = 0
            for dimlabel, idx in zip(data.dimlabels, range(len(data.dimlabels))):
                if dimlabel in stds:
                    if stds[dimlabel] == 0:
                        stds[dimlabel] = 0.0000001
                    dist += ((otherpoint[idx].chemshift - point[idx].chemshift) / stds[dimlabel]) ** 2
            normalized_euclidean_dist = sqrt(dist)
            if normalized_euclidean_dist <= chi2dist:
                region.append(otherpoint)
        return region

    def expand_cluster(self, point, neighbor_pts, cluster, data, stds):
        """Expand cluster by examining neighborhood region and add_member it to list of clusters.

        :param point: Represents single peak within peak list.
        :param list neighbor_pts: List of neighbor points.
        :param cluster: Represents particular cluster of points.
        :param data: Represents PeakList object.
        :type data: :class:`~ssc.physicalentities.PeakList`
        :type point: :class:`~ssc.physicalentities.Peak`
        :type cluster: :class:`~ssc.grouping.Cluster`
        :return: None
        :rtype: None
        """
        cluster.add_member(point)
        for n_pt in neighbor_pts:
            if n_pt not in self.visited:
                self.visited.append(n_pt)
                new_neighbor_pts = self.region_query(n_pt, data, stds)
                if len(new_neighbor_pts) >= self.minpts:
                    for nn_pt in new_neighbor_pts:
                        if nn_pt not in neighbor_pts:
                            neighbor_pts.append(nn_pt)

            for clstr in self.clusters.values():
                if n_pt not in clstr:
                    if n_pt not in cluster:
                        cluster.add_member(n_pt)

            # if p not yet member of any cluster, add_member p to cluster C
            if len(self.clusters) == 0:
                if n_pt not in cluster:
                    cluster.add_member(n_pt)

        self.clusters[self.count] = cluster
        self.global_clusters[self.datapath][self.count] = cluster

    def print_clusters(self):
        """Print clusters to the terminal.
        :return: None
        """
        # for cluster in self.clusters.values():
        #     print('Cluster #:', cluster.label)
        #     for p in cluster.members:
        #         print(p, p.assignment)
        #     print('===================================================')
        json_str = self._to_json()
        print(json_str)

        return json_str

    def write(self, filehandle, outformat="json"):
        """Write clustering results into file.

        :param filehandle: Results file pointer.
        :param str outformat: Rusults file format to use.
        :return: None
        :rtype: None
        """
        try:
            if outformat is "json":
                json_str = self._to_json()
                filehandle.write(json_str)
            else:
                filehandle.close()
                raise TypeError("Unknown file format.")
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')


    def _to_json(self):
        """Save list of clusters into JSON formatted string.

        :return:
        :rtype: str
        """
        clusters = []
        for cluster in self.global_clusters[self.datapath].values():
            cluster_dict = {"label": cluster.label,
                            "peaks": []}
            for idx, peak in enumerate(cluster.members):
                peak_dict = {"dimensions": peak.chemshiftslist,
                             "assignment": peak.assignmentslist,
                             "index": idx}
                cluster_dict["peaks"].append(peak_dict)
            clusters.append(cluster_dict)

        json_str = json.dumps(clusters, sort_keys=False, indent=4)
        return json_str
