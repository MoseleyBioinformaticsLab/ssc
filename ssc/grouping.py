#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.grouping
~~~~~~~~~~~~~~

This module implements algorithm for grouping peaks into proto-spin systems
using modified DBSCAN clustering algorithm.

"""

import sys
import numpy as np
import json
from scipy.stats import chi2
from math import sqrt

import bokeh.plotting as bkp
import bokeh.models as bkm

# Replace by this: http://stackoverflow.com/questions/13998901/generating-a-random-hex-color-in-python
COLOR_HEX = ['#FFFF00', '#FFC0CB', '#00CED1', '#9400D3', '#FA8072', '#20B2AA', '#00FFFF', '#1E90FF', '#FF8C00',
             '#800000', '#FF0000', '#B0E0E6', '#FFFFFF', '#708090', '#6495ED', '#696969', '#AFEEEE', '#FFE4B5',
             '#F5DEB3', '#6B8E23', '#FFFAFA', '#66CDAA', '#ADFF2F', '#EE82EE', '#00FF7F', '#A9A9A9', '#F0FFFF',
             '#FF69B4', '#FFEFD5', '#FFFACD', '#F0F8FF', '#D3D3D3', '#FFD700', '#FF7F50', '#FFDAB9', '#7B68EE',
             '#40E0D0', '#F5F5DC', '#FFA07A', '#FAA460', '#9932CC', '#696969', '#B22222', '#4B0082', '#FFFFF0',
             '#DA70D6', '#FFFFE0', '#778899', '#FDF5E6', '#6A5ACD', '#FFB6C1', '#00FFFF', '#DAA520', '#228B22',
             '#E9967A', '#7CFC00', '#BC8F8F', '#00FF00', '#FF1493', '#FF6347', '#FFE4E1', '#008000', '#DC143C',
             '#DCDCDC', '#0000FF', '#000000', '#FFA500', '#5F9EA0', '#B8860B', '#87CEFA', '#D3D3D3', '#D8BFD8',
             '#A52A2A', '#F5FFFA', '#00FA9A', '#FAEBD7', '#FFE4C4', '#FAF0E6', '#FFFAF0', '#191970', '#A0522D',
             '#2F4F4F', '#9ACD32', '#FF00FF', '#708090', '#8B0000', '#C71585', '#2F4F4F', '#FFF0F5', '#CD5C5C',
             '#808080', '#FFF8DC', '#90EE90', '#808080', '#BA55D3', '#00BFFF', '#B0C4DE', '#D2691E', '#F0FFF0',
             '#F0E68C', '#8FBC8F', '#98FB98', '#DDA0DD', '#EEE8AA', '#0000CD', '#32CD32', '#F8F8FF', '#D2B48C',
             '#FAFAD2', '#87CEEB', '#FF4500', '#800080', '#CD853F', '#F08080', '#DEB887', '#778899', '#FFF5EE',
             '#556B2F', '#808000', '#E6E6FA', '#3CB371', '#7FFF00', '#A9A9A9', '#DB7093', '#FFEBCD', '#7FFFD4',
             '#F5F5F5', '#00008B', '#8B008B', '#8A2BE2', '#BDB76B', '#C0C0C0', '#4169E1', '#483D8B', '#000080',
             '#E0FFFF', '#4682B4', '#2E8B57', '#006400', '#FF00FF', '#008080', '#ADD8E6', '#8B4513', '#FFDEAD',
             '#48D1CC', '#9370DB', '#008B8B']


class Cluster(object):
    """A cluster container for spin systems, i.e. groups of peaks (points)
    that are clustered together and belong to the same spin system.
    """

    def __init__(self, label, stds):
        """Initialize Cluster.

        :param str label: Cluster label (id).
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

    def get_members(self):
        """Return list of cluster members.

        :return: List of cluster members.
        :rtype: list
        """
        return self.members

    def has_member(self, member):
        """Test if member is inside members_list.

        :param member: Represents single peak within PeakList object.
        :type member: sass.pe.PhysicalEntities.Peak
        :return: True or False
        :rtype: bool
        """
        return member in self.members

    def get_x(self):
        """Get list of x coordinates of a cluster.

        :return: List of x coordinates.
        :rtype: list
        """
        return self.x

    def get_y(self):
        """Get list of y coordinates of a cluster.

        :return: List of y coordinates.
        :rtype: list
        """
        return self.y

    def get_mean_x(self):
        """Get mean of x coordinates of a cluster.

        :return: Mean of x coordinates.
        :rtype: float
        """
        return np.mean(self.get_x())

    def get_mean_y(self):
        """Get mean of y coordinates of a cluster.

        :return: Mean of y coordinates.
        :rtype: float
        """
        return np.mean(self.get_y())

    def get_desc(self):
        """Get description of a cluster member.

        :return: Member description.
        :rtype: list
        """
        return self.desc

    @property
    def representative(self):
        # points = np.array([[dim.chemshift for dim in point if dim.label in self.stds] for point in self.members])
        points = np.array([[dim.chemshift for dim in point] for point in self.members])
        representative = np.mean(points, axis=0)
        return representative


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
        self.noise = Cluster(label="-1", stds={})

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

        :param point: Represents single peak within PeakList object.
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

            for c in self.clusters.values():
                if not c.has_member(n_pt):
                    if not cluster.has_member(n_pt):
                        cluster.add_member(n_pt)

            # if p not yet member of any cluster, add_member p to cluster C
            if len(self.clusters) == 0:
                if not cluster.has_member(n_pt):
                    cluster.add_member(n_pt)

        self.clusters[self.count] = cluster
        self.global_clusters[self.datapath][self.count] = cluster

    def get_noise(self):
        """Get cluster of points (Peaks) that were clustered as 'noise'.

        :return: list of points that were clustered as 'noise'.
        :rtype: list
        """
        return self.clusters[-1].members

    def print_clusters(self):
        """Print clusters to the terminal.
        :return: None
        """
        for cluster in self.clusters.values():
            print('Cluster #:', cluster.label)
            for p in cluster.members:
                print(p, p.assignment)
            print('===================================================')

    def write(self, filehandle, outformat="json"):
        """

        :param filehandle:
        :param str outformat:
        :return:
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
        """Save list clusters into JSON formatted string.

        :return:
        """
        clusters = []
        for cluster in self.global_clusters[self.datapath].values():
            cluster_dict = {"name": cluster.label,
                            "peaks": []}
            for idx, peak in enumerate(cluster.members):
                peak_dict = {"dimensions": peak.chemshiftslist,
                             "assignment": peak.assignmentslist,
                             "index": idx}
                cluster_dict["peaks"].append(peak_dict)
            clusters.append(cluster_dict)

        json_str = json.dumps(clusters, sort_keys=False, indent=4)
        return json_str


    # def visualizeClusters(self, filename, title=""):
    #     """Visualize cluster in 2D space."""
    #     TOOLS = 'pan,box_zoom,box_select,resize,reset,wheel_zoom'
    #     bkp.output_file(filename)
    #     plot = bkp.figure(tools=TOOLS, title=title)
    #
    #     for cluster, color in zip(self.clusters, COLOR_HEX):
    #         source = bkp.ColumnDataSource(
    #             data=dict(
    #                 x=cluster.get_x(),
    #                 y=cluster.get_y(),
    #                 assignment=cluster.get_desc()
    #             )
    #         )
    #
    #         for point in cluster.members_list:
    #             if cluster.label == -1:
    #                 g = bkm.CircleX(x='x', y='y', line_color='#6666ee', fill_color='black', fill_alpha=0.5, size=12)
    #             else:
    #                 g = bkm.Circle(x='x', y='y', line_color='#6666ee', fill_color=color, fill_alpha=0.5, size=12)
    #
    #             gr = plot.add_glyph(source_or_glyph=source, glyph=g)
    #             g_hover = bkm.HoverTool(renderers=[gr],
    #                                     tooltips=[('x', '@x'), ('y', '@y'), ('assignment', '@assignment')])
    #             plot.add_tools(g_hover)
    #
    #         plot.text(cluster.get_mean_x(), cluster.get_mean_y(), text=[cluster.number],
    #                   text_color='black', text_align='center', text_font_size='10pt')
    #
    #     # bkp.show(plot)
    #     bkp.save(plot)
    #
    #     # def visualizeClusters(self):
    #     #     """Method that visualizes clusters in 2D space."""
    #     #
    #     #     TOOLS = 'pan,box_zoom,box_select,resize,reset,wheel_zoom'
    #     #
    #     #     output_file("points.html")
    #     #     plot = figure(tools=TOOLS)
    #     #
    #     #     for cluster, color in zip(self.clusters, COLOR_HEX):
    #     #         x_coordinates = cluster.get_x()
    #     #         y_coordinates = cluster.get_y()
    #     #
    #     #         plot.scatter(x_coordinates, y_coordinates, marker='circle',
    #     #                      line_color='#6666ee', fill_color=color, fill_alpha=0.5, size=12)      # legend=cluster.label
    #     #         plot.text(cluster.getXmean(), cluster.getYmean(), text=[cluster.clusterNumber],
    #     #                   text_color='black', text_align='center', text_font_size='10pt')
    #     #
    #     #     plot.scatter(self.Noise.get_x(), self.Noise.get_y(), marker='circle_x', line_color='#6666ee',
    #     #                  fill_color='black', fill_alpha=0.5, size=12)                              # legend=self.Noise.label
    #     #
    #     #     show(plot)