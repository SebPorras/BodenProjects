###############################################################################
# Date: 7/12
# Author: Sebastian Porras   
# Aims: Trying to implement python version of the POG graph.  
# To do: 
# 1) Create a sym node class for annotating   
# 2) Create a POG graph class  
###############################################################################

import json 
import numpy as np

class Edge():
    
    def __init__(self, start=None, end=None, edgeType=None, recip=None, backward=None, forward=None, weight = None) -> None:
        
        self._start = start 
        self._end = end 
        self._edgeType = edgeType 
        self._recip = recip 
        self._backward = backward 
        self._forward = forward
        self._weight = weight

    def setStart(self, start):
        self._start = start 

    def setEnd(self, end):
        if end == []:
            self._end = "*"
        else: 
            self._end = end 
    
    def setEdgeType(self, type):
        self._edgeType = type

    def setRecip(self, recip):
        self._recip = recip 
    
    def setBackward(self, back):
        self._backward = back
    
    def setForward(self, forward):
        self._forward = forward
    
    def setWeight(self, weight):
        self._weight = weight 


class SymNode():

    def __init__(self, name=None, value=None, edges=[]) -> None:
        self._name = name 
        self._value = value 
        self._edges = edges 

    def setName(self, name): 
        self._name = name 

    def setValue(self, val):
        self._value = val

    def addEdge(self, edge: Edge):
        self._edges.append(edge)

    def getName(self):
        return self._name

    def getValue(self):
        return self._value
        
    def getEdges(self):
        return self._edges
    

class POGraph():

    def __init__(self,version=None,indices=None, nodes=None,start=None,end = None,
        size=None, terminated=None, directed=None, name=None, isAncestor=None) -> None:
        
        self._version = version
        self._indices = indices
        self._nodes = nodes
        self._start = start
        self._end = end 
        self._size = size
        self._terminated = terminated 
        self._directed = directed
        self._name = name
        self._isAncestor = isAncestor 
    
    def setName(self, name: (str)): 
        self._name = name


    def POGraphFromJSON(self, jpog: dict, isAncestor: bool):
        """
        Parameters: 
            jpog(dict): JSON format of a POG for a sequence
            
            isAncestor(bool)
        """
        
        self.setName(jpog["Name"])
        self._size = jpog["Size"]
        self._directed = jpog["Directed"]
        self._terminated = jpog["Terminated"]
        self._isAncestor = isAncestor
        self._version = jpog["GRASP_version"]
        self._end = jpog["Ends"]
        self._start = jpog["Starts"]

        i_array = jpog["Indices"]
        indices = [i for i in i_array]
      
        self._indices = np.array(indices)

        adj = jpog["Adjacent"]

        if len(indices) != len(adj):
            raise RuntimeError + "File in incorrect format"
        
        #Assumption: If sequence is an extant, then the only possible edges 
        # will be between adjacent nodes (including virtual start/end 
        # nodes). Therefore, ancestor POGs can be created simply from edge 
        # indices. This is to avoid creating edges multiple times.   
        
        node_vals = jpog["Nodes"]
        nodes = []
        
        #first, instantiate all nodes 
        #Assumes that you only record outgoing edges
        for i in range(len(i_array)):
       
            node = SymNode()
                    
            node.setName(indices[i])
            node.setValue(node_vals[i]["Value"])
            
            edge = Edge()
            edge.setStart(indices[i])
            
            #set [] to * to signify end of seq 
            try: 
                edge.setEnd(adj[i][0]) 
            except IndexError: 
                edge.setEnd("*") 
            
            node.addEdge(edge)  
            nodes.append(node)
        
        self._nodes = np.array(nodes)
        
        if isAncestor:
        
            edgeInd = jpog["Edgeindices"]
            edges = jpog["Edges"]

            for i in range(len(jpog["Edgeindices"])):
                
                edge_info = edges[i]

                ##Current implementation for virtual start node, add it as an edge to the real (i.e index 0) start node 
                # This is based on assumption that start node is always indexed first based on indices array 
                if edgeInd[i][0] == -1:
                    
                    cor_node = self._nodes[0]

                    edge = Edge(start=edgeInd[i][0], end=edgeInd[i][1], edgeType=jpog["Edgetype"], 
                    recip=edge_info["Recip"], backward=edge_info["Backward"], forward=edge_info["Forward"], weight=edge_info["Weight"])

                    cor_node.addEdge(edge)

                else:
        
                    #find index for node where start of edge begins(returns a tuple with an array: very confusing )
                    node_loc = np.where(self._indices == edgeInd[i][0])[0][0]
                    
                    cor_node = self._nodes[node_loc]
                    
                    edge = Edge(start=edgeInd[i][0], end=edgeInd[i][1], edgeType=jpog["Edgetype"], 
                    recip=edge_info["Recip"], backward=edge_info["Backward"], forward=edge_info["Forward"], weight=edge_info["Weight"])

                    #check if there is a duplicate edge between adjacent residues 
                    #if yes, replace with more informative bidirectional edge 
                    adjacent_edge = cor_node._edges[0]             

                    if edgeInd[i][0] == adjacent_edge._start and edgeInd[i][1] == adjacent_edge._end:
                        
                        cor_node._edges[0] = edge
                        
                    else:
                        cor_node.addEdge(edge)
                
                
# with open("./python_structures/ASR.json", "r") as file:
#     data = json.load(file)


# extant = data["Input"]["Extants"]
# ancestors = data["Ancestors"]

# #print(ancestors[0].keys())

# case = extant[0]
# g = POGraph()
# g.POGraphFromJSON(case, False)

# print(g._name)

