#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.visualize
~~~~~~~~~~~~~

This module provides code to visualize spin system clusters.
"""

import os
import json
import random
from statistics import mean

import bokeh.plotting as bkp
import bokeh.models as bkm


random_color = lambda: random.randint(0, 255)


def visualize_clusters(clusters_path, x_idx, y_idx, x_label, y_label, title, plot_output_path=None):
    """Visualize cluster in 2D space.

    :param str clusters_path: Path to grouping results JSON file.
    :param str plot_output_path: Path where to save interactive HTML plot.
    :param str title: Plot title.
    :param int x_idx: Index of x dimension within peak.
    :param int y_idx: Index of y dimension within peak.
    :param str x_label: X axis label.
    :param str y_label: Y axis label.
    :return: None
    :rtype: :py:obj:`None`
    """
    TOOLS = "pan,box_zoom,box_select,resize,reset,wheel_zoom,save"

    if plot_output_path is None:
        plot_output_path = clusters_path + ".html"
        bkp.output_file(plot_output_path)
    else:
        bkp.output_file(plot_output_path)

    plot = bkp.figure(tools=TOOLS, title=title, toolbar_location="below")
    plot.title.text_font_size = "15pt"
    plot.title.align = "center"
    plot.xaxis.axis_label, plot.yaxis.axis_label = "{}, ppm".format(x_label), "{}, ppm".format(y_label)
    plot.xaxis.axis_label_text_font_size = "15pt"
    plot.yaxis.axis_label_text_font_size = "15pt"

    with open(clusters_path, "r") as infile:
        clusters_list = json.load(infile)

    for cluster in clusters_list:
        color = "#%02X%02X%02X" % (random_color(), random_color(), random_color())

        cluster_x_coords = []
        cluster_y_coords = []
        cluster_assignments = []

        for peak in cluster["peaks"]:
            cluster_x_coords.append(peak["dimensions"][x_idx])
            cluster_y_coords.append(peak["dimensions"][y_idx])
            cluster_assignments.append(peak["assignment"])

        cluster_x_centroid = mean(cluster_x_coords)
        cluster_y_centroid = mean(cluster_y_coords)

        source = bkp.ColumnDataSource(data=dict(x=cluster_x_coords, y=cluster_y_coords, assignment=cluster_assignments))

        for _ in cluster["peaks"]:
            if cluster["label"] == -1:
                g = bkm.CircleX(x='x', y='y', line_color='#6666ee', fill_color='black', fill_alpha=0.5, size=5)
            else:
                g = bkm.Circle(x='x', y='y', line_color='#6666ee', fill_color=color, fill_alpha=0.5, size=10)

            gr = plot.add_glyph(source_or_glyph=source, glyph=g)
            g_hover = bkm.HoverTool(renderers=[gr], tooltips=[('x', '@x'), ('y', '@y'), ('assignment', '@assignment')])
            plot.add_tools(g_hover)

        plot.text(x=cluster_x_centroid, y=cluster_y_centroid, text=[cluster["label"]],
                  text_color='black', text_align='center', text_font_size='10pt')

        if cluster["label"] >= 0:
            x_tolerance = cluster["stds"][x_label] * 4
            y_tolerance = cluster["stds"][y_label] * 4
            plot.ellipse(x=cluster_x_centroid, y=cluster_y_centroid, width=x_tolerance, height=y_tolerance, fill_color=None)

    print("saving plot to:", os.path.abspath(plot_output_path))
    bkp.save(plot)
