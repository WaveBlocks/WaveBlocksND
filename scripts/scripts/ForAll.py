#!/usr/bin/env python
"""The WaveBlocks Project

This file contains code some simple code to call a given
python script for a bunch of simulation result files.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 R. Bourquin
@license: Modified BSD License
"""

import argparse
from glob import glob
import subprocess as sp

from WaveBlocksND.FileTools import get_result_dirs, get_results_file
from WaveBlocksND import GlobalDefaults


# Move files with these extensions back to the simulationpath
# NOTE: It is the responsability of the code in the
#       'scriptcode' script to produce filenames that
#       do not yield collisions. Otherwise the files
#       will be overwritten without any warning!

save_extensions = [GlobalDefaults.ext_resultdatafile]

try:
    from WaveBlocksND.Interface import GraphicsDefaults
    save_extensions.append(GraphicsDefaults.output_format)
except:
    save_extensions.append(".png")


def execute_for_all(resultspath, command):
    """Call a given python script with the simulation results data file as last
    command line argument. The script and the data file are specified by (relative)
    file system paths.

    :param resultspath: The path where to look for simulation data.
    :param command: The python script that gets called for all simulations.
    """
    for simulationpath in get_result_dirs(resultspath):
        print(" Executing code for datafile in {}".format(simulationpath))

        # The file(s) with the simulation data
        resfiles = get_results_file(simulationpath)

        # Handle case where multiple hdf5 files are found
        if type(resfiles) is list:
            print("  Warning: more than one results file found!")
        else:
            resfiles = [resfiles]

        # For each file call the given script
        for resfile in resfiles:
            sp.call(command + [resfile])

            # Move newly created files back to the simulation path.
            # NOTE: It is the responsability of the code in the
            #       'scriptcode' script to produce filenames that
            #       do not yield collisions. Otherwise the files
            #       will be overwritten without any warning!
            for ext in save_extensions:
                for afile in glob("*" + ext):
                    sp.call(["mv", afile, simulationpath])




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--results",
                        type = str,
                        help = "Path where the results are.",
                        nargs = "?",
                        default = GlobalDefaults.path_to_results)

    parser.add_argument("-c", "--command",
                        type = str,
                        help = "The commands to execute.",
                        nargs = argparse.REMAINDER)


    args = parser.parse_args()

    print("** Will execute the command: " + " ".join(args.command))
    print("** for all files in '" + str(args.results) + "'")

    execute_for_all(args.results, args.command)

    print("Done")
