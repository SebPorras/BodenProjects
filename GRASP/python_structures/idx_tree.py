###############################################################################
# Date: 7/12
# Author: Sebastian Porras   
# Aims: This is an implementation of the IDx tree. Can be initialised from a 
# JSON object and created.Currently based on output from a joint reconstruction 
# Has very basic structures that allow most variables to be accessed or set.
# To do: 
# 1) Make a method that can parse this information into a nwk file so that 
# it could be read by other tree viewers.  
# 2) Investigate about adding in Branchpoint objects that store particular 
# info about what parents are and how far away they are etc. 
#################################################################################

import json
import numpy as np
import os

#print(os.getcwd())


class BranchPoint:

    def __init__(self, id, parent, dist, children, isleaf) -> None:
        
        self._id = id 
        self._parent = parent
        self._dist = dist
        self._children = children
        self._isLeaf = isleaf 

    def __str__(self) -> str:
     
        return f"Name: {self._id} \nParent Index: {self._parent}\nDistance To Parent {self._dist}\nChildren IDs: {self._children}\nLeaf: {self._isLeaf}"

class IdxTree:

    def __init__(self, n=0, useddistances=False):
        """
        n (int): number of branch points 
        useddistances (bool): True if distances have been used?
        """
        
        self._nBranches = n 
        self._branchpoints = np.array([None for i in range(n)])
        self._parents = np.array([None for i in range(n)])
        self._children = dict()
        self._indices = dict()
        
        if useddistances:
            self._distances = np.array([None for i in range(n)])
        else:
            self._distances = None


    def __str__(self) -> str:
        return f"Number of branchpoints: {self.getNBranches()}\nParents: {self.getParents()}\nChildren: {self.getChildren()}\nIndices: {self.getIndexes()}\nDistances: {self.getDistances()}"

    def getNBranches(self) -> int: return self._nBranches

    def getBranchpoints(self) -> np.array: return self._branchpoints

    def getParents(self) -> np.array: return self._parents

    def getChildren(self) -> dict: return self._children

    def getIndexes(self) -> dict: return self._indices

    def getDistances(self) -> np.array: return self._distances


    def getParent(self, CIdx: int) -> int:
        """
        Parameters:
            CIdx (int): Index of the child 

        Returns:
            int: Index of that child's parent 
        """       

        return self._parents[CIdx]

    def getIndex(self, name: str):

        return self._indices[name]


    def setParent(self, BIdx: int, PIdx: int):
        """
        Parameters: 
            BIdx (int): Branchpoint Index

            PIdx (int): Index of parent Branchpoint
        """

        self._parents[BIdx] = PIdx

    def setBranchPoint(self, BIdx: int, bp: BranchPoint):
        """
        Parameters: 
            BIdx (int): Branchpoint Index

            bp (Branchpoint): container to store information
        """

        self._branchpoints[BIdx] = bp
    
    def setDistance(self, BIdx: int, dist: float):
        """
        Parameters: 
            BIdx (int): Branchpoint Index

            dist (int): branch length of branchpoint
        """

        self._distances[BIdx] = dist

    def setIndex(self, BIdx: int, Pname: str):
        """
        Parameters: 
            BIdx (int): Branchpoint Index

            Pname (str): Name of ancestor or extant
        """
        
        self._indices[Pname] = BIdx

    def setChildren(self, PIdx: int, Children: np.array):
        """
        Parameters: 
            PIdx (int): Parent index

            Children (np.array): Array containing index of children
        """

        self._children[PIdx] = Children


    def distToParent(self, CIdx: int) -> float:
        """
        Parameters:
            CIdx (int): Index of the child 

        Returns:
            float: Distance to parent  
        """  

        return self._distances[CIdx]

    def isLeaf(self, BIDx) -> bool:
        """
        Parameters: 
            BIdx (int): Branchpoint Index
        """

        if self._children[BIDx] is None:
            return True
        else:
            return False

    def constructNwk(self) -> str:
        
        pass
        """
        nwk = "("
        i = 0
        if self.getChildren()[i] is None:
            return nwk + ");"
        else:
            i += 1
            self.constructNwk()
        """



# def IdxTreeFromJSON(json_file: json) -> IdxTree:
#     """
#     Instantiates a Idx Tree object from a Json file
    
#     Parameters:
#         json_file (json): JSON file object 
    
#     Returns:
#         IdxTree: 
#     """
#     with open(json_file, "r") as file:
#         json_file = json.load(file)

#     try: 
#         jdists = json_file["Input"]["Tree"]["Distances"]
#     except KeyError:
#         jdists = None

#     branch_num = json_file["Input"]["Tree"]["Branchpoints"]

#     tree = IdxTree(branch_num, jdists != None)

#     jlabels = json_file["Input"]["Tree"]["Labels"]
#     jparents = json_file["Input"]["Tree"]["Parents"]
    
#     #iterate at each branch point 
#     for i in range(branch_num):
        
#         try:

#             #index by parent name and assign what their parent branch point index is 
#             tree.setParent(i, jparents[i]) #tree._parents[i] = jparents[i]
            
#             if jdists is not None:

#                 #using same index as parent, order distances 
#                 tree.setDistance(i, jdists[i]) #tree._distances[i] = jdists[i] 
            
#             #index by parent name and assign branch point
#             tree.setIndex(i, jlabels[i]) #tree._index[jlabels[i]] = i
        
#         except RuntimeError:
#             print("Invalid JSON format")
    
#     #Next step is to record children of each parent
#     for PIdx in range(branch_num):
        
#         curr_children = []
    
#         for CIdx in range(branch_num):
            
#             if (tree.getParent(CIdx) == PIdx):
#                 curr_children.append(CIdx)
        
#         if len(curr_children) == 0:
#             ch_array = None
#         else:
#             ch_array = np.array(curr_children)
       
#         tree.setChildren(PIdx, ch_array)

#     for BIdx in range(branch_num):

#         bp = BranchPoint(jlabels[BIdx], tree.getParent(BIdx), jdists[BIdx],
#         tree.getChildren()[BIdx], tree.isLeaf(BIdx))
#         tree.setBranchPoint(BIdx, bp) 
    
#     return tree 

# test_tree = IdxTreeFromJSON("./python_structures/ASR.json")

# print(test_tree)

