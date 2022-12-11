###############################################################################
# Date: 9/12
# Author: Sebastian Porras   
# Aims: The POG tree will have two components, 1) the Idx tree which encodes 
# the topology of the tree. 2) The POG graph that will be stored according 
# to the corresponding index on the IDX tree so that both structures 
# are in agreement with each other. 
#################################################################################

import json 
import numpy as np
from pog_graph import *
from idx_tree import * 

class POGTree(object):
    """Contains an IdxTree and a POGraph for each branch point.
    All graphs and branchpoints are assigned the same index for 
    continuity
    """

    def __init__(self, idxTree:IdxTree=None, POGraphs=None) -> None:
        """Constructs instance of POGTree

        Parameters:
            idxTree(IdxTree)
            POGraph(POGraph)
        """
        
        self._idxTree = idxTree
        if POGraphs is None:
            self._POGraphs = dict()

    def getIdxTree(self) -> IdxTree: return self._idxTree

    def getPOGraphs(self) -> dict: return self._POGraphs

    def IdxTreeFromJSON(self, json_file: json):
        """Instantiates a IdxTree object from a Json file
        
        Parameters:
            json_file (json): JSON file object 
        """
        with open(json_file, "r") as file:
            json_file = json.load(file)

        try: 
            jdists = json_file["Input"]["Tree"]["Distances"]
        except KeyError:
            jdists = None

        branch_num = json_file["Input"]["Tree"]["Branchpoints"]

        tree = IdxTree(branch_num, jdists != None)

        jlabels = json_file["Input"]["Tree"]["Labels"]
        jparents = json_file["Input"]["Tree"]["Parents"]
        
        #iterate at each branch point 
        for i in range(branch_num):
            
            try:

                #index by parent name and assign what their parent 
                # branch point index is 
                tree.setParent(i, jparents[i]) 
                
                if jdists is not None:

                    #using same index as parent, order distances 
                    tree.setDistance(i, jdists[i])
                
                #index by parent name and assign branch point
                tree.setIndex(i, jlabels[i]) 
            
            except RuntimeError:
                print("Invalid JSON format")
        
        #Next step is to record children of each parent
        for PIdx in range(branch_num):
            
            curr_children = []
        
            for CIdx in range(branch_num):
                
                if (tree.getParent(CIdx) == PIdx):
                    curr_children.append(CIdx)
            
            if len(curr_children) == 0:
                ch_array = None
            else:
                ch_array = np.array(curr_children)
        
            tree.setChildren(PIdx, ch_array)

        for BIdx in range(branch_num):

            bp = BranchPoint(id=jlabels[BIdx], parent=tree.getParent(BIdx), 
            dist=jdists[BIdx],children=tree.getChildren()[BIdx],
            isleaf=tree.isLeaf(BIdx))
            
            tree.setBranchPoint(BIdx, bp) 
        
        self._idxTree = tree

    def POGraphFromJSON(self, JSON: json):
        """Instantiates a POGraph object from a Json file. Requires an 
        instance of an IdxTree so that graphs can be assigned to the correct 
        index
        
        Parameters:
            json_file (json): JSON file object 
        """
            
        with open(JSON, "r") as file:
            data = json.load(file)

        extants = data["Input"]["Extants"]
        ancestors = data["Ancestors"]

        for e in extants:

            g = POGraph()
            
            g.POGraphFromJSON(e, isAncestor=False)
          
            idx = self._idxTree.getIndex(g._name)

            self._POGraphs[idx] = g
        
        for a in ancestors:
            
            g = POGraph()
            g.POGraphFromJSON(a, isAncestor=True)

            idx = self._idxTree.getIndex(g._name)

            self._POGraphs[idx] = g


if __name__ == "__main__":
    

    poggers = POGTree()

    #First create the tree topology 
    poggers.IdxTreeFromJSON("./python_structures/ASR_big.json")

    #once the tree is created, you can assign POGs to branchpoints 
    poggers.POGraphFromJSON("./python_structures/ASR_big.json")

    print("IdxTree")
    print()
    print(poggers.getIdxTree())


    print()
    print("Graphs")
    print()
    print(poggers.getPOGraphs())

    #Example 1: Retrieve POG graph for ancestor zero 
    print()
    print("Tree indices")
    print(poggers.getIdxTree().getIndices())

    #Lets look at ancestor '6' -> can see its index is 9
    N6_pog = poggers.getPOGraphs()[9]
    print()
    print("POG summary info")
    print(N6_pog)
    print()
    #We can look for nodes that have multiple edges
    for n in N6_pog.getNodes():
     
        if len(n.getEdges()) >= 2:
            print()
            print(f"Sequence index: {n.getName()}")
            print(f"Number of edges: {len(n._edges)}")
            for e in n.getEdges():
                print(f"Start: {e.getStart()}, End: {e.getEnd()}")

    