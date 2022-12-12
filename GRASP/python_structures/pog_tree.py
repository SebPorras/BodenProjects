###############################################################################
# Date: 9/12
# Author: Sebastian Porras
# Aims: The POG tree will have two components, 1) the Idx tree which encodes
# the topology of the tree. 2) The POG graph that will be stored according
# to the corresponding index on the IDX tree so that both structures
# are in agreement with each other.
#################################################################################

import json
from pog_graph import *
from idx_tree import *


class POGTree(object):
    """Contains an IdxTree and a POGraph for each branch point.
    All graphs and branchpoints are assigned the same index for 
    continuity
    """

    def __init__(self, idxTree: IdxTree, POGraphs: dict) -> None:
        """Constructs instance of POGTree

        Parameters:
            idxTree(IdxTree)
            POGraph(POGraph)
        """

        self._idxTree = idxTree
        self._graphs = POGraphs

    def getIdxTree(self) -> IdxTree: return self._idxTree

    def getPOGraphs(self) -> dict: return self._graphs


def POGTreeFromJSON(json_path: str) -> POGTree:
    """Instantiates a POGraph object from a Json file. Requires an 
    instance of an IdxTree so that graphs can be assigned to the correct 
    index

    Parameters:
        json_path (json): path to JSON file
    """

    with open(json_path, "r") as file:
        data = json.load(file)

    tree = IdxTreeFromJSON(data)

    extants = data["Input"]["Extants"]
    ancestors = data["Ancestors"]

    graphs = dict()

    for e in extants:

        g = POGraphFromJSON(e, isAncestor=False)

        idx = tree.getIndex(g._name)

        graphs[idx] = g

    for a in ancestors:

        g = POGraphFromJSON(a, isAncestor=True)

        idx = tree.getIndex(g._name)

        graphs[idx] = g

    return POGTree(idxTree=tree, POGraphs=graphs)


if __name__ == "__main__":

    poggers = POGTreeFromJSON("./python_structures/ASR_big.json")

    print("IdxTree")
    print()
    print(poggers.getIdxTree())

    print()
    print("Graphs")
    print()
    print(poggers.getPOGraphs())

    # # Example 1: Retrieve POG graph for ancestor zero
    # print()
    # print("Tree indices")
    # print(poggers.getIdxTree().getIndices())

    # # Lets look at ancestor '6' -> can see its index is 9
    N6_pog = poggers.getPOGraphs()[36]
    # print()
    # print("POG summary info")
    print(N6_pog)

    print(poggers.getIdxTree().getIndex("XP_012291909.1"))

    # print()
    # # We can look for nodes that have multiple edges
    # for n in N6_pog.getNodes():

    #     if len(n.getEdges()) >= 2:
    #         print()
    #         print(f"Sequence index: {n.getName()}")
    #         print(f"Number of edges: {len(n._edges)}")
    #         for e in n.getEdges():
    #             print(f"Start: {e.getStart()}, End: {e.getEnd()}")
