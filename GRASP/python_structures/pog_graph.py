###############################################################################
# Date: 7/12
# Author: Sebastian Porras   
# Aims: Trying to implement python version of the POG graph.  
# To do: 
# 1) Create a sym node class for annotating   
# 2) Create a POG graph class  
###############################################################################

import numpy as np

class Edge(object):
    """Creates instance of an edge between two positions in a sequence. 
    Currently only works for bidirectional edges.
    
    """
    
    def __init__(self, start:int=None, end:int=None, edgeType:str=None,
    recip:bool=None, backward:bool=None, forward:bool=None, weight:float=None) -> None:
        """Constructs instance of an edge. 

        Parameters: 
            start(int): position of beginning of edge  
            end(int): position of end of edge
            edgeType(str): Currently only supports bidirectional edge
            recip(bool): ASK ABOUT THIS 
            backward(bool): Direction of edge
            forward(bool): Direction of edge
            weight(float): Support of the edge
        """
        
        self._start = start 
        self._end = end 
        self._edgeType = edgeType 
        self._recip = recip 
        self._backward = backward 
        self._forward = forward
        self._weight = weight

    def __str__(self) -> str:
        
        return(f"Start: {self._start}\nEnd: {self._end}\nType: \
        {self._edgeType}\nWeight: {self._weight}\nBackward: \
        {self._backward}\n Forward: {self._forward}\nRecip: \
        {self._recip}")

    def setStart(self, start: int): self._start = start 

    def setEnd(self, end: int):

        #assigns * to the virtual end node 
        if end == []:
            self._end = "*"
        else: 
            self._end = end 
    
    def setEdgeType(self, type: str): self._edgeType = type

    def setRecip(self, recip: bool): self._recip = recip 
    
    def setBackward(self, back: bool): self._backward = back
    
    def setForward(self, forward: bool): self._forward = forward
    
    def setWeight(self, weight: float): self._weight = weight 

    def getStart(self): return self._start

    def getEnd(self): return self._end

class SymNode(object):
    """Only implemented for output from joint reconstruction.
    Stores the most likely character at a given sequence position 
    and all of the edges at this position. 
    """

    def __init__(self, name:int=None, value:str=None, edges=None) -> None:
        """Constructs instance of the object. 
        
        Parameters:
            name(int): index position in sequence 
            value(str): Most likely amino acid based on joint reconstruction
            edges(list): Contains all outgoing edges at this position  
        """
        self._name = name 
        self._value = value 
        if edges is None: 
            self._edges = []

    def __str__(self) -> str:
        return f"Name: {self._name}\n Value: {self._value}\n\# of edges: {len(self._edges)}"

    def setName(self, name): self._name = name 

    def setValue(self, val): self._value = val

    def getName(self): return self._name

    def getValue(self): return self._value
        
    def getEdges(self): return self._edges

    def addEdge(self, e: Edge): self.getEdges().append(e)

class POGraph(object):
    """Representation of a sequence as a partial order graph (POG).
    Each position is assigned a SymNode which contain Edges. 
    """

    def __init__(self, version:str=None,indices:np.array=None, nodes:np.array=None,
    start:int=None, end:int=None, size:int=None, terminated:bool=None,
    directed:bool=None, name:str=None, isAncestor:bool=None):
        """Constructs instance of POGraph. 

        Parameters:
            version(str): Version of GRASP used for inference 
            indices(np.array): Contains all sequence position indices 
            nodes(np.array): contains SymNode objects for each position 
            start(int): Start index of sequence 
            end(int): End index of sequence 
            size(int): Size of sequence 
            terminated(bool): If the sequence has an end??
            directed(bool): If graph has order?? 
            name(str): Sequence ID 
            isAncestor(bool): Identifier for extant or ancestor 
        """
        
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

    def __str__(self) -> str:
        return f"Sequence ID: {self._name}\nSize: {self._size}\nStart: {self._start}\nEnd: {self._end}"
    
    def setName(self, name:str): self._name = name

    def setSize(self, size:int): self._size = size

    def setDirection(self, dir:bool): self._directed = dir

    def setTermination(self, ter:bool): self._terminated = ter

    def setAncestor(self, a:bool): self._isAncestor = a

    def setVersion(self, ver:str): self._version = ver

    def setEnd(self, end:int): self._end = end

    def setStart(self, start:int): self._start = start

    def getNodes(self): return self._nodes

    def POGraphFromJSON(self, jpog: dict, isAncestor: bool):
        """Takes JSON format of a POG and transcribes this information 

        Parameters: 
            jpog(dict): JSON format of a POG for a sequence
            
            isAncestor(bool): Only ancestors have multiple edges requiring 
            different fields to be read
        """
        
        #fill basic fields 
        self.setName(jpog["Name"])
        self.setSize(jpog["Size"])
        self.setDirection(jpog["Directed"])
        self.setTermination(jpog["Terminated"])
        self.setVersion(jpog["GRASP_version"])
        self.setStart(jpog["Starts"])
        self.setEnd( jpog["Ends"])

        indices = jpog["Indices"]
      
        self._indices = np.array(indices)

        adj = jpog["Adjacent"]

        if len(indices) != len(adj):
            raise RuntimeError + "File in incorrect format"
        
        #Assumption: If sequence is an extant, then the only possible edges 
        # will be between adjacent nodes (including virtual start/end 
        # nodes).
        
        node_vals = jpog["Nodes"]
        nodes = []
        
        #first, instantiate all nodes 
        #Assumes that you only record outgoing edges
        for i in range(len(indices)):
            
            node = SymNode(name=indices[i], value=node_vals[i]["Value"])
                    
            edge = Edge(start=indices[i])
            
            #set [] to * to signify end of seq 
            try: 
                edge.setEnd(adj[i][0]) 
            except IndexError: 
                edge.setEnd("*") 
            
            node.addEdge(edge)
            nodes.append(node)
     
        
        #print(len(nodes[10]._edges))
        #print(nodes)
        self._nodes = np.array(nodes)
        print(len(self._nodes))
        
        if isAncestor:
        
            edgeInd = jpog["Edgeindices"]
            edges = jpog["Edges"]

            for i in range(len(jpog["Edgeindices"])):
                
                edge_info = edges[i]

                # Current implementation for virtual start node, add it as an edge 
                # to the real start node (i.e index 0). This is based on assumption 
                # that start node is always indexed first based on indices array 
                if edgeInd[i][0] == -1:
                    
                    cor_node = self._nodes[0]

                    edge = Edge(start=edgeInd[i][0], end=edgeInd[i][1], 
                    edgeType=jpog["Edgetype"], recip=edge_info["Recip"],
                    backward=edge_info["Backward"], forward=edge_info["Forward"],
                    weight=edge_info["Weight"])

                    cor_node.addEdge(edge)

                else:

                    #edges are stored in a random order so need to find correct 
                    # node 
        
                    #find index for node where start of edge begins(returns a
                    # tuple with an array: very confusing)
                    node_loc = np.where(self._indices == edgeInd[i][0])[0][0]
                    
                    cor_node = self._nodes[node_loc]
                    
                    edge = Edge(start=edgeInd[i][0], end=edgeInd[i][1],
                    edgeType=jpog["Edgetype"], recip=edge_info["Recip"],
                    backward=edge_info["Backward"], forward=edge_info["Forward"],
                    weight=edge_info["Weight"])

                    #check if there is a duplicate edge between adjacent residues 
                    #if yes, replace with more informative bidirectional edge 
                    adjacent_edge = cor_node.getEdges()[0]             

                    if edgeInd[i][0] == adjacent_edge.getStart() and \
                    edgeInd[i][1] == adjacent_edge.getEnd():
                        
                        cor_node._edges[0] = edge
                        
                    else:
                        cor_node.addEdge(edge)
                      
                
