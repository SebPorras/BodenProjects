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
import numpy.typing as npt

#print(os.getcwd())


class BranchPoint(object):
    """Represents a branchpoint on a phylogenetic tree. Can contain
    information about the parents or children of that point and how 
    long that branch point is. 

    Future use -> can associate annotations with this branch point 
    """

    def __init__(self, id: str, parent: int, dist: float,
    children: npt.ArrayLike, isleaf: bool) -> None:
        """ Constructs instance of a branchpoint.

        Parameters: 
            id(str): Sequence ID 
            parent(int): Index of parent branchpoint
            dist(float): Distance to parent
            children(np.array): all children of branchpoint 
            isLeaf(bool): Marker for extant sequence             
        """
        
        self._id = id 
        self._parent = parent
        self._dist = dist
        self._children = children
        self._isLeaf = isleaf 

    def __str__(self) -> str:
     
        return (f"Name: {self._id}\n\
        Parent Index: {self._parent}\n\
        Distance To Parent {self._dist}\n\
        Children IDs: {self._children}\n\
        Leaf: {self._isLeaf}")

    def getId(self) -> str: return self._id

    def getParent(self) -> int: return self._parent

    def getDist(self) -> float: return self._dist 

    def getChildren(self): return self._children

    def isLeaf(self) -> bool: return self._isLeaf 

class IdxTree(object):

    def __init__(self, n=0, useddistances=False):
        """
        Parameters: 
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
        return f"Number of branchpoints: {self.getNBranches()}\nParents: {self.getParents()}\nChildren: {self.getChildren()}\nIndices: {self.getIndices()}\nDistances: {self.getDistances()}"

    def getNBranches(self) -> int: return self._nBranches

    def getBranchpoints(self) -> np.array: return self._branchpoints

    def getParents(self) -> np.array: return self._parents

    def getChildren(self) -> dict: return self._children

    def getIndices(self) -> dict: return self._indices

    def getDistances(self) -> np.array: return self._distances

    def getParent(self, CIdx: int) -> int:
        """
        Parameters:
            CIdx (int): Index of the child 

        Returns:
            int: Index of that child's parent 
        """       

        return self._parents[CIdx]

    def getIndex(self, name: str) -> int:
        """
        Parameters:
            name (str): sequence name 

        Returns:
            int: index of branchpoint for that sequence
        
        """

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

    def setChildren(self, PIdx: int, Children: npt.ArrayLike):
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