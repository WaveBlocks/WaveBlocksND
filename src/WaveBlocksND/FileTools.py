"""The WaveBlocks Project

This file contains various functions for finding and retrieving
the files that contain parameter settings and simulation results.
Note: The terms 'path' and 'ID' are used as synonyms here. Each simulation
ID is just the basename of the path or the configuration file.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import os

import GlobalDefaults as GD


def get_result_dirs(path):
    """Lists all simulations (IDs) that can be found under the given path.

    :param path: The filesystem path under which we search for simulations.
    :return: A list of simulation IDs.
    """
    dirs = [ os.path.join(path, dir) for dir in os.listdir(path) ]
    return dirs


def get_parameters_file(path):
    """Search for a configuration file containing the simulation parameters under a given path.
    Note that in case there are more than one .py file under the given path we just return the
    first one found!

    :parameter path: The path under which we search for a configuration file.
    :return: The path (filename) of the configuration file.
    """
    parameters_file = None

    for file in os.listdir(path):
        if file[-3:] == ".py":
            parameters_file = file
            break

    if parameters_file is None:
        raise IOError("No configuration .py file found!")

    parameters_file = os.path.join(path, parameters_file)
    return parameters_file


def get_results_file(path):
    """Search for a file containing the simulation results under a given path.
    Note that in case there are more than one .hdf5 file under the given path
    we just return the first one found!

    :param path: The path under which we search for a output file.
    :return: The path (filename) of the output file.
    """
    results_file = None

    for file in os.listdir(path):
        # Should we allow extensions .hdf and .h5 too?
        if file[-5:] == ".hdf5":
            results_file = file
            break

    if results_file is None:
        raise IOError("No results .hdf5 file found!")

    results_file = os.path.join(path, results_file)
    return results_file


def get_number_simulations(path):
    """Get the number of simulations at hand below the given path.

    :param path: The path under which we search for a output file.
    :return: The number of simulations result directories.
    """
    ids = get_result_dirs(path)
    return len(ids)


def name_contains(name, pattern):
    """Checks if a simulation ID contains a given pattern.

    :param name: The full simulation ID.
    :param pattern: The pattern in question.
    :return: A boolean answer.
    """
    return pattern in name


def gather_all(stringlist, pattern):
    """Collects all simulation IDs which contain a specific pattern from a given list.

    :param stringlist: A list with the simulation IDs
    :param pattern: The pattern
    :return: A list of simulation IDs that contain the given pattern.
    """
    gathered = [ s for s in stringlist if name_contains(s, pattern) ]
    return gathered


def compare_by(namea, nameb, pattern, ldel=GD.kvp_ldel, mdel=GD.kvp_mdel, rdel=GD.kvp_rdel, as_string=True):
    """Compare two simulation IDs with respect to a (numerical) value in the ID.

    :param namea: The first name in the comparison
    :param nameb: The second name in the comparison
    :param pattern: The pattern whose (numerical) value is used for sorting
    :param ldel: Left delimiter of the pattern
    :param mdel: Middle delimiter of the pattern
    :param rdel: Right delimiter of the pattern
    :param as_string: Determines if the values for ``pattern`` get converted to floats
    :return: A boolean answer if the IDs are the same w.r.t the pattern.
    """
    part1 = namea.partition(ldel + pattern + mdel)
    part2 = nameb.partition(ldel + pattern + mdel)

    part1 = part1[0:2] + part1[2].partition(rdel)
    part2 = part2[0:2] + part2[2].partition(rdel)

    # Convert values to float for comparisons
    if not as_string:
        part1 = part1[0:2] + (float(part1[2]),) + part1[3:]
        part2 = part2[0:2] + (float(part2[2]),) + part2[3:]

    return part1[2] == part2[2]


def group_by(stringlist, pattern, ldel=GD.kvp_ldel, mdel=GD.kvp_mdel, rdel=GD.kvp_rdel, as_string=True):
    """Groups simulation IDs with respect to a pattern.

    :param stringlist: A list with the simulation IDs
    :param pattern: The pattern used for grouping
    :param ldel: Left delimiter of the pattern
    :param mdel: Middle delimiter of the pattern
    :param rdel: Right delimiter of the pattern
    :param as_string: Determines if the values for ``pattern`` get converted to floats. Not used here.
    :return: A list of groups of simulation IDs.
    """
    tmp = [ s.partition(ldel + pattern + mdel) for s in stringlist ]
    tmp = [ s[0:2] + s[2].partition(rdel) for s in tmp ]

    distinct_vals = set([ s[2] for s in tmp ])

    groups = [ [] for i in xrange(len(distinct_vals)) ]

    for item in tmp:
        for index, val in enumerate(distinct_vals):
            # Yes, we compare the strings here and not the floats
            # to avoid representation errors when converting to floats.
            if item[2] == val:
                # Concatenate the fragments again
                groups[index].append( reduce(lambda x,y: x+y, item) )
                break

    return groups


def intersect_by(lista, listb, pattern, ldel=GD.kvp_ldel, mdel=GD.kvp_mdel, rdel=GD.kvp_rdel, as_string=True):
    """Find the intersection of two lists containing simulation IDs.

    :param lista: A first list with the simulation IDs
    :param listb: A second list with the simulation IDs
    :param pattern: The pattern whose numerical value is used for sorting
    :param ldel: Left delimiter of the pattern
    :param mdel: Middle delimiter of the pattern
    :param rdel: Right delimiter of the pattern
    :param as_string: Determines if the values for ``pattern`` get converted to floats
    :return: A sorted list of simulation IDs.
    """
    result = []

    for namea in lista:
        for nameb in listb:
            if compare_by(namea, nameb, pattern, ldel=ldel, mdel=mdel, rdel=rdel, as_string=as_string):
                result.append(namea)
                result.append(nameb)

    # Remove possible duplicates
    result = list(set(result))
    return result


def sort_by(stringlist, pattern, ldel=GD.kvp_ldel, mdel=GD.kvp_mdel, rdel=GD.kvp_rdel, as_string=False):
    """Sorts simulation IDs with respect to a (numerical) value in the ID.

    :param stringlist: A list with the simulation IDs
    :param pattern: The pattern whose (numerical) value is used for sorting
    :param ldel: Left delimiter of the pattern
    :param mdel: Middle delimiter of the pattern
    :param rdel: Right delimiter of the pattern
    :param as_string: Determines if the values for ``pattern`` get converted to floats
    :return: A sorted list of simulation IDs.
    """
    tmp = [ s.partition(ldel + pattern + mdel) for s in stringlist ]
    tmp = [ s[0:2] + s[2].partition(rdel) for s in tmp ]

    if not as_string:
        # Convert to float and append numeric values to the splitted IDs
        tmp = [ s + (float(s[2]),) for s in tmp ]
    else:
        # Use string in comparison, allows sorting
        tmp = [ s + (s[2],) for s in tmp ]

    # Define a costum order relation
    def compare(x,y):
        if x[-1] > y[-1]:
            return 1
        elif y[-1] > x[-1]:
            return -1
        else:
            return 0

    # Sort w.r.t. the numerical value
    tmp.sort(cmp=compare)

    # Remove numeric value and concatenate the fragments again
    f = lambda x,y: x+y
    sorted_list = [ reduce(f, s[:-1]) for s in tmp ]

    return sorted_list


def get_max_by(stringlist, pattern, ldel=GD.kvp_ldel, mdel=GD.kvp_mdel, rdel=GD.kvp_rdel, as_string=False):
    """Get the maximum of a list with simulation IDs with respect to a (numerical) value in the ID.
    This is just a simple convenience function so that the user needs not to remember if the
    sort order is ascending or descending which plays no role for iteration.

    :param stringlist: A list with the simulation IDs
    :param pattern: The pattern whose (numerical) value is used for sorting
    :param ldel: Left delimiter of the pattern
    :param mdel: Middle delimiter of the pattern
    :param rdel: Right delimiter of the pattern
    :param as_string: Determines if the values for ``pattern`` get converted to floats
    :return: A sorted list of simulation IDs.
    """
    sortedlist = sort_by(stringlist, pattern, ldel=ldel, mdel=mdel, rdel=rdel, as_string=as_string)
    return sortedlist[-1]


def get_min_by(stringlist, pattern, ldel=GD.kvp_ldel, mdel=GD.kvp_mdel, rdel=GD.kvp_rdel, as_string=False):
    """Get the minimum of a list with simulation IDs with respect to a (numerical) value in the ID.
    This is just a simple convenience function so that the user needs not to remember if the
    sort order is ascending or descending which plays no role for iteration.

    :param stringlist: A list with the simulation IDs
    :param pattern: The pattern whose (numerical) value is used for sorting
    :param ldel: Left delimiter of the pattern
    :param mdel: Middle delimiter of the pattern
    :param rdel: Right delimiter of the pattern
    :param as_string: Determines if the values for ``pattern`` get converted to floats
    :return: A sorted list of simulation IDs.
    """
    sortedlist = sort_by(stringlist, pattern, ldel=ldel, mdel=mdel, rdel=rdel, as_string=as_string)
    return sortedlist[0]
