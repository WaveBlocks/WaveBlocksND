#!/usr/bin/env python
"""The WaveBlocks Project

Simple script to run several simulations with a given set of parameter files.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os
import shlex
import logging
import concurrent.futures
import subprocess


def list_configurations(path):
    r"""Search for configuration files containing the simulation parameters
    under a given path and return a list of all.

    :parameter path: The path under which we search for configuration files.
    """
    configurations = []

    for afile in os.listdir(path):
        afile = os.path.join(path, afile)
        if os.path.isfile(afile) and afile.endswith('.py'):
            configurations.append(afile)

    return configurations


class Command:
    r"""
    """

    def __init__(self, commandstring):
        self._command = self._get_lambda(shlex.split(commandstring))
        self._commandstring = commandstring

    def _get_lambda(self, commandstring):
        def lambda_c(configfile, outputpath):
            U = []
            for c in commandstring:
                if c == "CONFIGFILE":
                    U.append(configfile)
                elif c == "OUTPUTPATH":
                    U.append(outputpath)
                else:
                    U.append(c)
            return U
        return lambda_c

    def __repr__(self):
        return self._commandstring

    @property
    def command(self):
        return self._command



class Job:
    r"""
    """

    def __init__(self, configfile, resultspath, commands):
        self._configfile = configfile
        self._resultspath = resultspath
        self._commands = commands

    def __repr__(self):
        return repr(self._commands)

    @property
    def resultspath(self):
        return self._resultspath

    @property
    def configfile(self):
        return self._configfile

    def commands(self, outputpath):
        configfile = self._configfile
        return [ command.command(configfile, outputpath) for command in self._commands ]



def run_job(job):
    r"""Run the given job object.
    """
    print("Running simulation: {}".format(job.configfile))

    configfile = job.configfile
    configname = os.path.basename(configfile).split('.py')[0]

    outputpath = os.path.join(job.resultspath, configname)

    stdoutlogfile = os.path.join(outputpath, "simulation_"+str(configname)+"_stdout.log")
    stderrlogfile = os.path.join(outputpath, "simulation_"+str(configname)+"_stderr.log")

    # Make sure the ouput directory exists
    os.makedirs(outputpath, mode=0o700, exist_ok=True)

    # Open log files and execute commands
    with open(stdoutlogfile, 'a') as stdout:
        with open(stderrlogfile, 'a') as stderr:
            for command in job.commands(outputpath):
                subprocess.call(command,
                                shell=False,
                                stdout=stdout,
                                stderr=stderr)
    return None



def batch_loop(jobs, max_workers=4):
    r"""Run jobs in parallel.

    :param jobs: List of runnable jobs objects.
    :param max_workers: Maximal number of jobs run in parallel.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:

        F = [ executor.submit(run_job, job) for job in jobs ]

        for future in concurrent.futures.as_completed(F):
            try:
                future.result()
            except Exception:
                logging.exception('Code generated an exception')
                continue





if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--configurations",
                        type = str,
                        required = True,
                        help = "Path to the 'configuration' directory.")

    parser.add_argument("-r", "--resultspath",
                        type = str,
                        help = "Path to the 'results' directory.",
                        default = '.')

    parser.add_argument("-m", "--maxworkers",
                        type = int,
                        help = "Maximal number of parallel jobs.",
                        default = 2)

    args = parser.parse_args()


    commands = [
        "cp CONFIGFILE OUTPUTPATH",
        "Main.py CONFIGFILE -r OUTPUTPATH",
        "ComputeNorms.py -r OUTPUTPATH",
        "ComputeEnergies.py -r OUTPUTPATH",
        "ComputeAutocorrelations.py -r OUTPUTPATH",
        "PlotNorms.py -r OUTPUTPATH",
        "PlotEnergies.py -r OUTPUTPATH",
        "PlotAutocorrelations.py -r OUTPUTPATH"
    ]


    # Paths
    fpath = os.path.abspath(args.configurations)
    rpath = os.path.abspath(args.resultspath)

    if not (os.path.exists(fpath) and os.path.isdir(fpath)):
        raise ValueError("Invalid configuration path: " + str(fpath))

    if not (os.path.exists(rpath) and os.path.isdir(rpath)):
        raise ValueError("Invalid results path: " + str(rpath))

    # Read off commands
    C = [ Command(c) for c in commands ]

    # List all configuration files
    F = list_configurations(fpath)

    # Set up jobs
    J = [ Job(f, rpath, C) for f in F ]

    # Batch run
    print("Running {} simulations from: {}".format(len(F), fpath))
    print("Putting results into: {}\n".format(rpath))
    print("-------------------")
    batch_loop(J, max_workers=args.maxworkers)
