#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Routines for grouping peaks from single peak list into spin system clusters.

This package includes the following modules:

``ssc``
    This module provides the :class:`~ssc.ssc.SpinSystemCreator` class which
    manages creation of spin system clusters.

``registration``
    This module provides the :func:`~ssc.registration.run_registration` function
    that executes registration algorithm executable and reports results of the
    registration algorithm as python :py:class:`dict`.

``grouping``
    This module provides the :class:`~ssc.grouping.Cluster` class for storing peaks
    that belong to the same spin system and :class:`~ssc.ssc.DBSCAN` that groups
    peaks into spin system clusters.

``physicalentities``
    This module provides interfaces to represents physical entities such as
    :class:`~ssc.physicalentities.Dimension`, :class:`~ssc.physicalentities.Peak`,
    :class:`~ssc.physicalentities.PeakList`.

``peaklistparsers``
    This module provides different parsers to create :class:`~ssc.physicalentities.PeakList`
    from peak lists of different formats.

``visualize``
    This module provides code to visualize spin system clusters as interactive graph.
"""

__version__ = "0.1.0"
