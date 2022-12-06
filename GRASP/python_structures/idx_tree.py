import json
import numpy as np


def readJSON(json_file):
    with open(json_file, "r") as file:
        data = json.load(file)

    return data

class Idx_Tree:

    def __init__(self, n: int, useddistances: bool):
        """
        n (int): number of branch points 
        useddistances (bool): True if distances have been used?
        """
        
        self.branchpoints = np.array([None for i in range(n)])
        self.parents = np.array([None for i in range(n)])
        #have used dict here instead of matrix due new int[branch_len][] data type being difficult in python 
        self.children = dict()
        self.index = dict()
        
        if useddistances:
            self.distances = np.array([None for i in range(n)])
        else:
            self.distances = None

def fromJSON(json_file):

    try: 
        jdists = json_file["Input"]["Tree"]["Distances"]
    except KeyError:
        jdists = None

    branch_num = json_file["Input"]["Tree"]["Branchpoints"]

    tree = Idx_Tree(branch_num, jdists != None)

    jlabels = json_file["Input"]["Tree"]["Labels"]
    jparents = json_file["Input"]["Tree"]["Parents"]
    
    #iterate at each branch point 
    for i in range(branch_num):
        
        try:
            #index by parent name and assign what their parent branch point index is 
            tree.parents[i] = jparents[i]
            
            if jdists is not None:
                #using same index as parent, order distances 
                tree.distances[i] = jdists[i] 
            #index by parent name and assign branch point 
            tree.index[jlabels[i]] = i
        
        except RuntimeError:
            print("Invalid JSON format")
    
    #Next step is to record children of each parent
    for parent in range(branch_num):

        curr_children = []
        
        for child in range(branch_num):
            
            if tree.parents[child] == parent:
                curr_children.append(child)
        
        if len(curr_children) == 0:
            continue
        else:
            #need to check if it's alright to use dictionary 
            #also unsure if instances of parents without children is recorded in a particular way 
            tree.children[parent] = np.array([None for i in range(len(curr_children))]) 
            for i in range(len(curr_children)):
                tree.children[parent][i] = curr_children[i]

    return tree 

data = readJSON("ASR.json")
print(data["Input"]["Extants"][0].keys())
test_tree = fromJSON(data)
print(test_tree.children)





