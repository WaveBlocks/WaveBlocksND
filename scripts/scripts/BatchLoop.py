#!/usr/bin/env python
"""The WaveBlocks Project

Simple script to run several simulations with a given set of parameter files.

@author: R. Bourquin
@copyright: Copyright (C) 2016 R. Bourquin
@license: Modified BSD License
"""

import argparse
import os

import concurrent.futures
import subprocess


class job:

    def __init__(self, configfile, resultspath):
        self._configfile = configfile
        self._resultspath = resultspath

    def resultspath(self):
        return self._resultspath

    def configfile(self):
        return self._configfile

    def commands(self, resultspath):
        configfile = self._configfile
        return [
            ['Main.py', str(configfile), '-r', resultspath],
            ['ComputeNorms.py', '-r', resultspath],
            ['ComputeEnergies.py', '-r', resultspath],
            ['ComputeAutocorrelations.py', '-r', resultspath],
            ['PlotNorms.py', '-r', resultspath],
            ['PlotEnergies.py', '-r', resultspath],
            ['PlotAutocorrelations.py', '-r', resultspath]
        ]



def list_configurations(path):
    configurations = [
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_1D_p.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_1D_f.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_1D_p_ih.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_1D_p_nsd.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_1D_p_stationary_groundstate.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_2D_p.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_2D_p_stationary_groundstate.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_ring_2D_p.py',
        '/u/raoulb/WaveBlocksND/examples/harmonic_oscillators/harmonic_tube_2D_p.py'
    ]
    return configurations



def run_job(job):
    print("Running simulation: " + job.configfile())

    resultspath = job.resultspath()
    configfile = job.configfile()
    configname = os.path.basename(configfile).split('.py')[0]

    with open("call_"+str(configname)+"_stdout.log", 'w') as stdout:
        with open("call_"+str(configname)+"_sdterr.log", 'w') as stderr:

            outputpath = os.path.join(resultspath, configname)
            os.mkdir(outputpath)

            for command in job.commands(outputpath):
                subprocess.call(command,
                                shell=False,
                                stdout=stdout,
                                stderr=stderr)

    return True



def batch_loop(jobs, max_workers=4):

    with concurrent.futures.ThreadPoolExecutor(max_workers = max_workers) as executor:

        F = {executor.submit(run_job, job) : job for job in jobs}

        for future in concurrent.futures.as_completed(F):
            try:
                future.result()
            except Exception:
                print('Code generated an exception')
                continue
            else:
                print('Code terminated properly')
                pass




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-j", "--job",
                        type = str,
                        help = "The job configuration file")

    parser.add_argument("-c", "--configurations",
                        type = str,
                        help = "Path to the 'configuration' directory")

    parser.add_argument("-r", "--results",
                        type = str,
                        help = "Path to the 'results' directory.",
                        default = '.')

    parser.add_argument("-m", "--maxworkers",
                        type = int,
                        help = "Maximal number of parallel jobs.",
                        default = 4)

    args = parser.parse_args()


    # Set up jobs
    J = []
    for c in list_configurations('.'):
        J.append(job(c, './results'))


    # Batch run
    batch_loop(J, max_workers=4)
