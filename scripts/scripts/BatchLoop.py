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
        afile = os.path.abspath(os.path.join(path, afile))
        if os.path.isfile(afile) and afile.endswith('.py'):
            configurations.append(afile)

    return configurations



def parse_commands(commands):
    def get_lambda(command):
        def lambda_c(configfile, outputpath):
            U = []
            for c in command:
                if c == "CONFIGFILE":
                    U.append(configfile)
                elif c == "OUTPUTPATH":
                    U.append(outputpath)
                else:
                    U.append(c)
            return U
        return lambda_c

    L = [get_lambda(shlex.split(command)) for command in commands ]
    return L



class Job:
    r"""
    """
    def __init__(self, configfile, resultspath, commands):
        self._configfile = configfile
        self._resultspath = resultspath
        self._commands = commands

    @property
    def resultspath(self):
        return self._resultspath

    @property
    def configfile(self):
        return self._configfile

    def commands(self, outputpath):
        configfile = self._configfile
        return [ command(configfile, outputpath) for command in self._commands ]



def run_job(job):
    r"""Run the given job object.
    """
    print("Running simulation: " + job.configfile)

    resultspath = job.resultspath
    configfile = job.configfile
    configname = os.path.basename(configfile).split('.py')[0]

    with open("simulation_"+str(configname)+"_stdout.log", 'w') as stdout:
        with open("simulation_"+str(configname)+"_stderr.log", 'w') as stderr:

            outputpath = os.path.join(resultspath, configname)

            for command in job.commands(outputpath):
                subprocess.call(command,
                                shell=False,
                                stdout=stdout,
                                stderr=stderr)
    return True



def batch_loop(jobs, max_workers=4):
    r"""Run jobs in parallel.

    :param jobs: List of runnable jobs objects.
    :param max_workers: Maximal number of jobs run in parallel.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:

        F = {executor.submit(run_job, job) : job for job in jobs}

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
                        help = "Path to the 'configuration' directory")

    parser.add_argument("-r", "--resultspath",
                        type = str,
                        help = "Path to the 'results' directory.",
                        default = './results')

    parser.add_argument("-m", "--maxworkers",
                        type = int,
                        help = "Maximal number of parallel jobs.",
                        default = 2)

    args = parser.parse_args()


    commands = [
        "mkdir -p OUTPUTPATH",
        "Main.py CONFIGFILE -r OUTPUTPATH",
        "ComputeNorms.py -r OUTPUTPATH",
        "ComputeEnergies.py -r OUTPUTPATH",
        "ComputeAutocorrelations.py -r OUTPUTPATH",
        "PlotNorms.py -r OUTPUTPATH",
        "PlotEnergies.py -r OUTPUTPATH",
        "PlotAutocorrelations.py -r OUTPUTPATH"
    ]


    # Read off commands
    C = parse_commands(commands)

    # List all configuration files
    F = list_configurations(args.configurations)

    # Set up jobs
    J = [ Job(f, args.resultspath, C) for f in F ]

    # Batch run
    batch_loop(J, max_workers=args.maxworkers)
