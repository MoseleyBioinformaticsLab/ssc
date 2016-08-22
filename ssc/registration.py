#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ssc.registration
~~~~~~~~~~~~~~~~~~

This module implements helper function that runs registration algorithm (C++ binary).

"""

import subprocess
import sys
import re
import json


def run_registration(execpath, inputlistpath, rootlistpath, regresultpath, *args):
    """Run registration algorithm executable and save results.

    :param str execpath: Path to registration algorithm executable.
    :param str inputlistpath: Path to input peak list.
    :param str rootlistpath: Path to root peak list.
    :param str regresultpath: Path where save registration results.
    :param args: Extra command-line arguments necessary to run registration executable.
    :return: Results dict.
    :rtype: dict
    """
    results = {}
    numberpattern = re.compile("[+-]?\d+(\.\d+)?")
    cmdarguments = [execpath, inputlistpath, rootlistpath] + list(args)
    print(cmdarguments)
    output = subprocess.check_output(cmdarguments)
    lines = output.decode("utf-8").splitlines()

    if not lines:
        return results
    else:
        with open(regresultpath, "r") as infile:
            results = json.load(infile)

    return results

if __name__ == "__main__":

    script = sys.argv.pop(0)
    ipl = sys.argv.pop(0)
    rpl = sys.argv.pop(0)
    resultspath = sys.argv.pop(0)

    r = run_registration('./bin/calculate_registration', ipl, rpl, resultspath, '-autoassign', '-noi', '-save',
                         resultspath, '-dim', 'HN', 'N', ':', 'HN', 'N')
