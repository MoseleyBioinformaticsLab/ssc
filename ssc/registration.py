#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.registration
~~~~~~~~~~~~~~~~~~

This module implements helper function that runs the 
registration algorithm.

"""

import subprocess
import json


def run_registration(regalg_executable, input_peaklist, root_peaklist, *args):
    """Run registration algorithm executable and save results.

    :param str regalg_executable: Path to the registration algorithm executable.
    :param str input_peaklist: Path to input peak list.
    :param str root_peaklist: Path to root peak list.
    :param args: Command-line arguments necessary to run registration executable.
    :return: Dictionary containing registration results.
    :rtype: :py:class:`dict`
    """
    results = {}
    command = "{} {} {} ".format(regalg_executable, input_peaklist, root_peaklist) + " ".join(list(args))
    output = subprocess.check_output(command, shell=True)
    lines = output.decode("utf-8")

    if not lines:
        return results
    else:
        try:
            results = json.loads(lines)
        except ValueError:
            return results
    return results
