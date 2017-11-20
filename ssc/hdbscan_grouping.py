#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.hdbscan_grouping
~~~~~~~~~~~~~~~~~~~~

This module implements algorithm for grouping peaks into proto-spin systems
using modified DBSCAN clustering algorithm.
"""

import sys

import hdbscan
import numpy as np
import bokeh.plotting as bkp
import bokeh.models as bkm

from . import peaklistparsers as plp
from . import grouping


def run_hdbscan(peaklist_path, result_path, spectrum_type, dim_labels, root_dim_labels, pl_format, min_num_pts=2):

    dbc = grouping.DBSCAN(peaklist_path, min_num_pts)

    peaklist = plp.parse(peaklist_path, spectrum_type, dim_labels, pl_format)
    clusterer = hdbscan.HDBSCAN(min_num_pts)

    cluster_labels = list(clusterer.fit_predict(peaklist.peaklistdf[root_dim_labels]))
    clusters = {}

    for peak, label in zip(peaklist, cluster_labels):
        if label >= 0:
            label += 1
        clusters.setdefault(int(label), [])
        clusters[int(label)].append(peak)

    for label in clusters.keys():
        clstr_obj = grouping.Cluster(label=label, stds={})
        for peak in clusters[label]:
            clstr_obj.add_member(peak)
        dbc.global_clusters[peaklist_path][label] = clstr_obj

    with open(result_path, "w") as outfile:
        dbc.write(outfile)

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


def visualize_clusters(peaklistpath, plotfilepath, title, x, y):
    """Visualize cluster in 2D space."""

    TOOLS = "pan,box_zoom,box_select,resize,reset,wheel_zoom,save"
    bkp.output_file(plotfilepath)
    plot = bkp.figure(tools=TOOLS, title=title)

    plot.xaxis.axis_label, plot.yaxis.axis_label = "{}, ppm".format(x), "{}, ppm".format(y)

    for cluster, color in zip(grouping.DBSCAN.global_clusters[peaklistpath].values(), COLOR_HEX):
        x_coords = cluster.coordinates(x)
        y_coords = cluster.coordinates(y)
        assignments = cluster.assignments()
        source = bkp.ColumnDataSource(data=dict(x=x_coords, y=y_coords, assignment=assignments))

        for point in cluster.members:
            if cluster.label == -1:
                g = bkm.CircleX(x='x', y='y', line_color='#6666ee', fill_color='black', fill_alpha=0.5, size=5)
            else:
                g = bkm.Circle(x='x', y='y', line_color='#6666ee', fill_color=color, fill_alpha=0.5, size=10)

            gr = plot.add_glyph(source_or_glyph=source, glyph=g)
            g_hover = bkm.HoverTool(renderers=[gr],
                                    tooltips=[('x', '@x'), ('y', '@y'), ('assignment', '@assignment')])
            plot.add_tools(g_hover)

        plot.text(np.mean(x_coords), np.mean(y_coords), text=[cluster.label],
                  text_color='black', text_align='center', text_font_size='10pt')
    # bkp.show(plot)
    bkp.save(plot)

if __name__ == "__main__":
    # python3 -m ssc.hdbscan_grouping datasets/sparky/HNcoCACB/jr19_hncocacb.txt jr19_hncocacb.json HNcoCACB H,N,CA/CB-1 H,N sparky

    script = sys.argv.pop(0)
    peaklistfilepath = sys.argv.pop(0)
    resultfilepath = sys.argv.pop(0)
    spectrumtype = sys.argv.pop(0)
    dimlabels = sys.argv.pop(0).split(",")
    rootdimlabels = sys.argv.pop(0).split(",")
    plformat = sys.argv.pop(0)

    run_hdbscan(peaklist_path=peaklistfilepath,
                result_path=resultfilepath,
                spectrum_type=spectrumtype,
                dim_labels=dimlabels,
                root_dim_labels=rootdimlabels,
                pl_format=plformat)

    visualize_clusters(peaklistfilepath, "{}.html".format(resultfilepath), spectrumtype, *rootdimlabels)
