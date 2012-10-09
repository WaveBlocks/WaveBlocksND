"""The WaveBlocks Project

This file contains code some simple code to call a given
python script for a bunch of simulation result files.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import sys
from glob import glob
import subprocess as sp
import time

from WaveBlocksND.FileTools import get_result_dirs, get_results_file
from WaveBlocksND import GlobalDefaults

try:
    from GraphicsDefaults import output_format
except:
    output_format = ".png"


def execute_for_all(resultspath, scriptcode, scriptargs):
    """Call a given python script with the simulation results data file as first
    command line argument. The script and the data file are specified by (relative)
    file system paths.

    :param resultspath: The path where to look for simulation data.
    :param scriptcode: The python script that gets called for all simulations.
    :param scriptargs: Optional (constant) arguments to the script.
    """
    # UNBUFFERED timelog file
    timelog = open("logtime_forall", "w", 0)

    for simulationpath in get_result_dirs(resultspath):
        starttime = time.time()
        timelog.writelines(["New script starting at: " + time.ctime(starttime) + "\n",
                            "Executing code for datafile in " + simulationpath + "\n"])
        print(" Executing code for datafile in " + simulationpath)

        # The file with the simulation data
        afile = get_results_file(simulationpath)

        # Call the given script
        sp.call(["python", scriptcode, afile] + scriptargs)

        endtime = time.time()
        timelog.writelines(["Script now finished at timestamp: " + time.ctime(endtime) + "\n",
                            "Script took " + str(endtime-starttime) + " second to run\n\n"])

        # Move plots away if any
        for afile in glob("*"+output_format):
            sp.call(["mv", afile, simulationpath])

    timelog.close()




if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Call syntax is:")
        print("---------------")
        print(" ForAll.py script.py")
        print(" ForAll.py resultpath script.py")
        print(" ForAll.py resultpath script.py script_args")
        sys.exit()
    elif len(sys.argv) == 2:
        # python ForAll.py script.py
        scriptcode = sys.argv[1]
        resultspath = GlobalDefaults.path_to_results
        scriptargs = []
    elif len(sys.argv) == 3:
        # python ForAll.py resultpath script.py
        scriptcode = sys.argv[2]
        resultspath = sys.argv[1]
        scriptargs = []
    elif len(sys.argv) > 3:
        # python ForAll.py resultpath script.py script_args
        scriptcode = sys.argv[2]
        resultspath = sys.argv[1]
        scriptargs = sys.argv[3:]

    print("Will execute the code from '" + scriptcode + "' for all files in '" + resultspath + "'")
    print("with optional arguments " + str(scriptargs))

    execute_for_all(resultspath, scriptcode, scriptargs)

    print("Done")
