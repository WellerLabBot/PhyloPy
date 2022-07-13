"""

    Phylogenetic Analysis Program
    Written By: Eliza Diggins

"""
### IMPORTS ###
import matplotlib.pyplot as plt
import logging as log
import os
import numpy as np
import pandas as pd
import json
from matplotlib.collections import LineCollection as LC

### CORE VARS ###
_supported_tree_filetypes = [".dnd"]


### CLASSES ###
class Tree():
    def __init__(self, raw, name):
        """
        Builds the tree given the raw input
        :param raw: string, either filename or a raw input of the data in tuples.
        :type raw: str or tuple
        """
        ### Intro Logging ###
        log.debug("BioPython:Phylogeny:Tree:__init__:DEBUG: Creating Tree(%s)" % str(raw))

        ### Sanitizing input ###
        if isinstance(raw, str):
            # This should be a filetype
            log.debug(
                "BioPython:Phylogeny:Tree:__init__:DEBUG: 'raw' was type (str), attempting to read file at %s." % raw)
            if os.path.isfile(raw) and raw[-4:] in _supported_tree_filetypes:
                # This is a valid tree
                log.debug(
                    "BioPython:Phylogeny:Tree:__init__:DEBUG: Found file %s with extension %s. Attempting read." % (
                    raw, raw[:-4]))
                self.tree = read_tree(raw)
            else:
                # This is not a valid tree file
                log.warning("BioPython:Phylogeny:Tree:__init__:WARNING: Input of type str was not a file, assuming trivial tree.")
                self.tree = raw
        elif isinstance(raw, tuple):
            # This should be a valid raw input
            log.debug(
                "BioPython:Phylogeny:Tree:__init__:DEBUG: 'raw' was type (tuple), attempting to construct from raw.")
            self.tree = raw  # Building the tree variable
        else:
            log.error(
                "BioPython:Phylogeny:Tree:__init__:ERROR: Failed to recognize input 'raw' of type %s. CONSTRUCTION FAILED." % type(
                    raw))
            del self

        ### Generating Data ###
        self.leaves = self.read_leaves()
        self.name = name

    def read_leaves(self):
        # Reads the leaves of the input
        return read_leaves(self.tree)
    
    def find_clades(self,distance=0.01):
        """
        Finds all of the clades at the given distance from the origin.
        :param distance: 
        :return: 
        """
        self.clades = recur_clades(self.tree,target_distance=distance)
        return self.clades

### Sub Functions ###
def find_continents_from_country(countries:list,dataset_location:str)->list:
    """
    Converts all of the countries into their respective continents.
    :param countries: The list of countries to search for.
    :param dataset_location: The location of the translation JSON.
    :return: List of the continents
    """
    ### Loading the file ###
    try:
        dataset = json.load(open(dataset_location))
        dataset = pd.read_json(json.dumps(dataset,ensure_ascii=False).encode("latin-1").decode("cp1252"))
    except Exception:
        log.error("BioPython:Phylogeny:find_continents_from_country:ERROR: Failed to open %s as json."%dataset_location)
        return [np.nan for country in countries]

    output_list = []

    for country in countries:
        try:
            output_list.append(dataset.loc[dataset["country"]==country,"continent"].item())
        except:
            output_list.append("")

    return output_list

### CORE FUNCTIONS ###
def recur_clades(tree_tup,recursion_number=0,distance=0.0,target_distance=1,n_clades=0):
    """
    Finds clades from the tuple by the following algorithm:
        For each branch in the phylogeny, check the distance. If distance > distance ---> This is already a clade ----> We create a tree and add it to the clades.
        
        Else: Explore that branch recursively until we find a clade.

    :param tree_tup: The tuple form of the tree to analyze.
    :param recursion_number: The recursion number, always start from 0.
    :param distance: The distance of the origin of the branch.
    :param target_distance: The distance we want to target for the clade formation.
    :return: Return the list of tree objects each representing the given clade.
    """
    ### INTRO LOGGING ###
    log.debug("BioPython:Phylogeny:recur_clades:DEBUG: Computing clades of %s at distance %s. Recursion Number = %s."%(tree_tup,distance,recursion_number))

    ### VARS ###
    CLADES = []
    ### CHECKING IF THE BRANCH FORMS A CLADE ###
    for element in tree_tup:
        # recursively checking each element in the tree tuple at this level.

        ## Grabbing the element's distance
        try:
            element_distance = float(element[-1]) + distance
        except ValueError:
            log.warning("BioPython:Phylogeny:recur_clades:WARNING: Failed to find element distance for element %s."%str(element))
            element_distance = None

        if element_distance:
            # The element distance exists.
            if element_distance >= target_distance:
                # We have identified a branch or element that is longer than the target length and is therefore a clade.
                log.debug("BioPython:Phylogeny:recur_clades:DEBUG: %s was found to be clade number %s."%(str(element[0]),n_clades+len(CLADES)+1))
                ## Turning the element into a valid clade ##
                if isinstance(element[0],str):
                    ## The element is a string.
                    CLADES.append(Tree(element[0],name="Clade_%s"%(n_clades+len(CLADES)+1))) # adds all of the clades

                elif isinstance(element[0],tuple):
                    ### The element is a new branch ###
                    CLADES.append(Tree(element[0],name="Clade_%s"%(n_clades+len(CLADES)+1)))
            else:
                # We have a branch that needs to be further analyzed.
                if isinstance(element[0],str):
                    # The element is a string and it must never get that far, so we pass it.
                    pass
                else:
                    CLADES += recur_clades(element[0],recursion_number=recursion_number+1,distance=element_distance,target_distance=target_distance,n_clades=n_clades+len(CLADES))

    return CLADES
    
def read_tree(filepath: str) -> tuple:
    """
    Reads a .dnd file and returns a tree
    :param filepath: The filepath to open.
    :return: The tree in the correct format as a tuple.
    """
    # Intro Logging #
    log.debug('BioPython:Phylogeny:read_tree:DEBUG: Reading tree from %s.' % filepath)

    # Grabbing the data #
    with open(filepath, "r+") as file:
        data = file.read()

    return format_list(read_string(data))


def read_string(string, recursion_number=0) -> list:
    """
    Reads the string and converts it to the correct tuple format. This is done as follows:

    1. Remove redundant parenthesis and split on ":" if necessary
    :param string: The input string.
    :return: Tuple of the correct format.
    """
    # Intro Logging #
    log.debug("BioPython:Phylogeny:read_string:DEBUG: Recursion: %s; string: %s." % (recursion_number, string))

    if recursion_number == 0:
        # We need to manage the header.
        temp_string = string[1:-2]  # Removing the ( .... ); from the outsides of the file.
        working_string = temp_string  # setting the working string.

        # spliting on delimiter
        lst = find_correct_comma(working_string)  # break by the correct ","

        for index, element in enumerate(lst):
            if element[0] == "(":  # we have more work to do
                lst[index] = read_string(element, recursion_number=recursion_number + 1)
            else:
                lst[index] = element.split(":")

        return lst
    else:
        # we need to manage normal body.
        temp_string = string.rsplit(":", 1)  # Splits on the last occurrence of the colon.

        # Managing parenthesis
        if temp_string[0][0] + temp_string[0][-1] == "()":
            temp_string[0] = temp_string[0][1:-1]

        temp_string[0] = find_correct_comma(temp_string[0])

        for index, element in enumerate(temp_string[0]):
            if element[0] == "(":  # we have more work to do
                temp_string[0][index] = read_string(element, recursion_number=recursion_number + 1)
            else:
                temp_string[0][index] = element.split(":")

        return temp_string


def format_list(inp_list: list, recursion_number: int = 0) -> tuple:
    """
    Re-formats the tree list into a tuple form.
    :param inp_list: The input list to work with.
    :return: Returns the tuple form of the tree.
    """
    # Intro logging
    log.debug("BioPython:Phylogeny:format_list:DEBUG: Formatting list %s. (recursion number=%s)." % (
    inp_list, recursion_number))

    for index, element in enumerate(inp_list):
        if isinstance(element, list):  # the constituent elements are lists
            if isinstance(element[0], str):  # The first element is a string ---> this is a base level.
                inp_list[index] = tuple(element)
            else:
                inp_list[index] = tuple(format_list(element, recursion_number=recursion_number + 1))

    return tuple(inp_list)

def read_leaves(tree:tuple,recursion_number=0)->list:
    """
    Input a dictionary of the form {{{leaf_1,leaf_2},leaf_3},{leaf_4,leaf_5}} for example, and returns leaves in a printable order.
    :param recursion_number: Simply keeps track of the number of recursions occurring.
    :param tree: The dictionary to parse. Leaves must be strings or ints.
    :return: list of leaf names.
    """
    log.debug("BioPython:Phylogeny:read_leaves:DEBUG: Recursion: %s; Tree: %s."%(recursion_number,tree))

    ### Checking for string input ###
    if isinstance(tree,str):
        # The input is a single trivial tree:
        return [tree]

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

if __name__ == '__main__':
    log.basicConfig(level=log.DEBUG)
    #t = Tree("/home/ediggins/BioInformatics/BioPython/tree.dnd",name="test")
    t = phylo.Tree("/media/Mercury/SRR_SALMONELLA_FILES/salmonella.dnd",name="test")
    recur_clades(t.tree)