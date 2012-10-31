"""The WaveBlocks Project

This file contains some tree manipulating functions.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

def validate_path(path):
    r"""Validate the path and resolve upward references.
    """
    parts = path.split("/")
    parts = filter(lambda p: len(p)>0, parts)
    P = []
    for part in parts:
        if part == ".":
            continue
        if part == "..":
            try:
                P.pop()
                continue
            except IndexError:
                raise ValueError("Invalid path: "+path)
        P.append(part)
    newpath = "/" + reduce(lambda a,b: a+"/"+b, P)
    return newpath


def split_path(path):
    r"""Split the path at each '/' and return the parts.
    """
    path = validate_path(path)
    parts = path.split("/")
    parts = filter(lambda x:len(x)>0, parts)
    return parts


def is_reference(value):
    r"""Check if a value is a reference.
    """
    if type(value) is str and value.startswith("@"):
        return True


def resolve_reference_path(path, value):
    r"""Resolve the references into one single path.
    """
    value = value[1:]
    if value.startswith("/"):
        # Case of absolute path
        return value
    elif value.startswith(".."):
        # Case of relative path
        return path +"/"+ value
    else:
        raise ValueError("Invalid relative reference: "+value)


def get_node(tree, path):
    r"""Retrieve a node from the tree given its path.
    """
    parts = split_path(path)
    node = tree
    for part in parts:
        node = node[part]
    return node


def get_value(tree, path):
    r"""Retrieve the value from the tree given the patch of its node.
    """
    p = path
    v = get_node(tree, p)
    while is_reference(v):
        p = resolve_reference_path(p, v)
        v = get_node(tree, p)
    return v


def merge_tree(treea, path, treeb):
    r"""Merge two trees. Plug one tree into the other below a given path.
    """
    parts = split_path(path)

    node = treea
    for part in parts:
        if not node.has_key(part):
            node[part] = {}

        node = node[part]
        print(node.keys())

    for k,v in treeb.iteritems():
        node[k] = v

    return treea
