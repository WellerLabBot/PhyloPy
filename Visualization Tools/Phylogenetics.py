"""

    Visualization Tools For Phylogenetic Analysis
                Written By: Eliza Diggins

"""
### IMPORTS ###
import matplotlib.pyplot as plt
import logging as log
import os
import numpy as np
from matplotlib.collections import LineCollection as LC
### CORE VARS ###
_supported_tree_filetypes = [".dnd"]

### CLASSES ###
class Tree():
    def __init__(self,raw,name):
        """
        Builds the tree given the raw input
        :param raw: string, either filename or a raw input of the data in tuples.
        :type raw: str or tuple
        """
        ### Intro Logging ###
        log.debug("BioPython:Phylogenetics:Tree:__init__:DEBUG: Creating Tree(%s)"%raw)

        ### Sanitizing input ###
        if isinstance(raw,str):
            # This should be a filetype
            log.debug("BioPython:Phylogenetics:Tree:__init__:DEBUG: 'raw' was type (str), attempting to read file at %s."%raw)
            if os.path.isfile(raw) and raw[-4:] in _supported_tree_filetypes:
                # This is a valid tree
                log.debug("BioPython:Phylogenetics:Tree:__init__:DEBUG: Found file %s with extension %s. Attempting read."%(raw,raw[:-4]))
                self.tree = read_tree(raw)
            else:
                # This is not a valid tree file
                log.error("BioPython:Phylogenetics:Tree:__init__:ERROR: File %s was not recognized as either a file or a valid extension. CONSTRUCTION FAILED.")
                del self
        elif isinstance(raw,tuple):
            # This should be a valid raw input
            log.debug("BioPython:Phylogenetics:Tree:__init__:DEBUG: 'raw' was type (tuple), attempting to construct from raw.")
            self.tree = format_list(read_string(raw))# Building the tree variable
        else:
            log.error("BioPython:Phylogenetics:Tree:__init__:ERROR: Failed to recognize input 'raw' of type %s. CONSTRUCTION FAILED."%type(raw))
            del self

        ### Generating Data ###
        self.leaves = self.read_leaves()
        self.name = name


    def read_leaves(self):
        # Reads the leaves of the input
        return read_leaves(self.tree)

    def plot_tree(self,save=True,labels=True,text_offset = 0.01):
        """
        Plots the phylogenetic tree.
        :param save: True to save, False to show
        :param labels: True to show labels on the leaves.
        :return: None
        """
        # Intro logging
        log.debug('BioPython:Phylogenetics:Tree:plot_tree:DEBUG: Plotting tree of %s')

        # Building the figure
        fig = plt.figure()
        ax1 = fig.add_subplot(111)


        # Computing the line segments
        lines = compute_segments(self.tree)

        ax1.add_collection(LC(lines[0]+lines[1]))

        if labels:
            # We are actually going to add labels
            for index,label in enumerate(self.leaves):
                ax1.text(float(lines[2][index][0])*(1+text_offset),lines[2][index][1],label)

        ax1.autoscale()
        plt.show()



### CORE FUNCTIONS ###
def read_tree(filepath:str)->tuple:
    """
    Reads a .dnd file and returns a tree
    :param filepath: The filepath to open.
    :return: The tree in the correct format as a tuple.
    """
    # Intro Logging #
    log.debug('BioPython:Phylogenetics:read_tree:DEBUG: Reading tree from %s.'%filepath)

    # Grabbing the data #
    with open(filepath,"r+") as file:
        data = file.read()

    return format_list(read_string(data))

def read_string(string,recursion_number=0) -> list:
    """
    Reads the string and converts it to the correct tuple format. This is done as follows:

    1. Remove redundant parenthesis and split on ":" if necessary
    :param string: The input string.
    :return: Tuple of the correct format.
    """
    # Intro Logging #
    log.debug("BioPython:Phylogenetics:read_string:DEBUG: Recursion: %s; string: %s."%(recursion_number,string))


    if recursion_number == 0:
        # We need to manage the header.
        temp_string = string[1:-2] # Removing the ( .... ); from the outsides of the file.
        print(temp_string)
        working_string = temp_string # setting the working string.
        
        # spliting on delimiter
        lst = find_correct_comma(working_string) # break by the correct ","
        
        for index,element in enumerate(lst):
            if element[0] == "(": # we have more work to do
                lst[index] = read_string(element,recursion_number=recursion_number+1)
            else:
                lst[index] = element.split(":")

        return lst
    else:
        # we need to manage normal body.
        temp_string = string.rsplit(":",1) # Splits on the last occurrence of the colon.

        # Managing parenthesis
        if temp_string[0][0]+temp_string[0][-1] == "()":
            temp_string[0] = temp_string[0][1:-1]
        
        temp_string[0] = find_correct_comma(temp_string[0])

        for index,element in enumerate(temp_string[0]):
            if element[0] == "(": # we have more work to do
                temp_string[0][index] = read_string(element,recursion_number=recursion_number+1)
            else:
                temp_string[0][index] = element.split(":")

        return temp_string

def format_list(inp_list:list,recursion_number:int = 0)->tuple:
    """
    Re-formats the tree list into a tuple form.
    :param inp_list: The input list to work with.
    :return: Returns the tuple form of the tree.
    """
    # Intro logging
    log.debug("BioPython:Phylogenetics:format_list:DEBUG: Formatting list %s. (recursion number=%s)."%(inp_list,recursion_number))


    for index,element in enumerate(inp_list):
       if isinstance(element,list): # the constituent elements are lists
           if isinstance(element[0],str): # The first element is a string ---> this is a base level.
               inp_list[index] = tuple(element)
           else:
               inp_list[index] = tuple(format_list(element,recursion_number=recursion_number+1))

    return tuple(inp_list)





def read_leaves(tree:tuple,recursion_number=0)->list:
    """
    Input a dictionary of the form {{{leaf_1,leaf_2},leaf_3},{leaf_4,leaf_5}} for example, and returns leaves in a printable order.
    :param recursion_number: Simply keeps track of the number of recursions occurring.
    :param tree: The dictionary to parse. Leaves must be strings or ints.
    :return: list of leaf names.
    """
    log.debug("BioPython:Phylogenetics:read_leaves:DEBUG: Recursion: %s; Tree: %s."%(recursion_number,tree))
    output_list = [] # Create blank storage list.
    for element in tree:
        if isinstance(element[0],(str)): # is the element a string or an int?
            output_list.append(element[0])
        elif isinstance(element[0],(tuple)): # the branch forms another sub tree,
            output_list += read_leaves(element[0],recursion_number=recursion_number+1)

    return output_list

def find_correct_comma(string):
    """
    Finds the comma that isn't wrapped by parenthesis
    :param string: The string to split.
    :return: split case.
    """
    array = np.array([char for char in string])
    comma_locations = np.where(array==",")[0]

    splits = []
    for comma_loc in comma_locations:
        if len(np.where(array[:comma_loc]=="(")[0]) == len(np.where(array[:comma_loc]==")")[0]):
            splits.append(comma_loc)

    return [string[:splits[0]]]+[string[splits[i]+1:splits[i+1]] for i in range(0,len(splits)-1)] + [string[splits[-1]+1:]]


### COMPUTATION FUNCTIONS ###


def compute_segments(tup,origin=(0,0),recursion_number = 0,w_unit=1,v_unit=1):
    """
    Computes the vertical and horizontal lines for the given tree
    :param tup: The tree object to pass through the function.
    :param origin: The origin from which to begin computations
    :param recursion_number: The recursion number
    :return: Lines
    """
    log.debug("BioPython:Phylogenetics:compute_segments:DEBUG: Computing segments for %s. Recursion: %s."%(tup,recursion_number))


    # Creating data holders
    vertical_lines = [] # vertical lines are added as (line). Line has format [(x_0,y_0),(x_1,y_1)]
    horizontal_lines = [] # The horizontal lines are added.
    branch_endpoints = [] # Markers for the branch end points.

    # Computing leaves remaining
    leaf_counts = []
    for element in tup:
        if isinstance(element[0],str): # This is a single leaf
            leaf_counts.append(1)
        else: # This is not a single leaf, so we compute the number of remaining leaves.
            leaf_counts.append(len(read_leaves(element[0])))


    # Creating the vertical line
    v_line_length = sum(leaf_counts) # This gives the height of the vertical line.
    vertical_lines.append([(origin[0],origin[1]+(v_line_length/2) - 0.5*(leaf_counts[-1])),(origin[0],origin[1]-(v_line_length/2)+0.5*(leaf_counts[0]))])

    # Computing offsets
    v_offsets = [origin[1]-(v_line_length/2)+ sum(leaf_counts[:j]) + 0.5*leaf_counts[j] for j in range(len(leaf_counts))]

    for index,element in enumerate(tup): # We now cycle through again
        branch_length = float(element[1])
        horizontal_lines.append([(origin[0],v_offsets[index]),
                                 (origin[0]+branch_length,v_offsets[index])])

        if not isinstance(element[0],str):
            new_lines = compute_segments(element[0],origin=(origin[0]+branch_length,v_offsets[index]),recursion_number=recursion_number+1)
            vertical_lines += new_lines[0]
            horizontal_lines += new_lines[1]
            branch_endpoints += new_lines[2]
        else:
            branch_endpoints.append((origin[0]+branch_length,v_offsets[index]))

    ### Managing unit changes
    vertical_lines = [[(u[0][0]*w_unit,u[0][1]*v_unit),(u[1][0]*w_unit,u[1][1]*v_unit)] for u in vertical_lines]
    horizontal_lines = [[(u[0][0] * w_unit, u[0][1] * v_unit), (u[1][0] * w_unit, u[1][1] * v_unit)] for u in
                      horizontal_lines]
    branch_endpoints = [(u[0]*w_unit,u[1]*v_unit) for u in branch_endpoints]

    return [vertical_lines,horizontal_lines,branch_endpoints]






### VISUALIZATION FUNCTIONS ###



if __name__ == '__main__':
    #log.basicConfig(level=log.DEBUG)
    t = Tree("/home/ediggins/BioInformatics/BioPython/tree.dnd",name="test")
    t.plot_tree()