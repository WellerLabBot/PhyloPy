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
from Modules import Phylogeny as phylo

### CORE VARS ###
_supported_tree_filetypes = [".dnd"]

### SUB PROCESSES ###
def get_random_colors(n: int, name: str = "hsv") -> np.array:
    """
    Returns an array of HSV colors evenly distributed over n values.
    :param n: The number of values to use in the array
    :param name: The name of the colormap to use
    :return: Array containing colors.
    """
    return [plt.cm.get_cmap(name, n)(i) for i in range(n)]



### CORE FUNCTIONS ###

def read_leaves(tree:tuple,recursion_number=0)->list:
    """
    Input a dictionary of the form {{{leaf_1,leaf_2},leaf_3},{leaf_4,leaf_5}} for example, and returns leaves in a printable order.
    :param recursion_number: Simply keeps track of the number of recursions occurring.
    :param tree: The dictionary to parse. Leaves must be strings or ints.
    :return: list of leaf names.
    """
    log.debug("BioPython:Phylogenetics:read_leaves:DEBUG: Recursion: %s; Tree: %s."%(recursion_number,tree))
    if isinstance(tree,str):
        return [tree]
    output_list = [] # Create blank storage list.
    for element in tree:
        if isinstance(element[0],(str)): # is the element a string or an int?
            output_list.append(element[0])
        elif isinstance(element[0],(tuple)): # the branch forms another sub tree,
            output_list += read_leaves(element[0],recursion_number=recursion_number+1)

    return output_list


def compute_segments(tup,origin=(0,0),recursion_number = 0,w_unit=1,v_unit=1,clades=None,clade_colors=None):
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
    vcolors = []
    hcolors = []

    # grabbing leaves
    tup_leaves = read_leaves(tup)

    # Computing clades list
    if clades:
        """
        We have been given clades, so we will use clade coloring. This works as follows:
        
        If all of the leaves in this tuple are a part of the clade, then the whole set of lines is colored the clade_coloring.
        
        If not, we check each future branch and color it if possible. Finally, each other branch is colored black.
        
        """
        clade_leaves = [clade.leaves for clade in clades]
        if any(all(leaf in clade_leaf for leaf in tup_leaves) for clade_leaf in clade_leaves):
            # There is a tuple wide match
            index = [tup_leaves[0] in clade_leaf for clade_leaf in clade_leaves].index(True)

            tuple_color = clade_colors[index]
        else:
            tuple_color = None

        # Now we check if the tuple color failed. If it did, we need to check element by element.
        if not tuple_color:
            branch_colors = []
            for element in tup:
                element_leaves = read_leaves(element[0])
                if any(all(leaf in clade_leaf for leaf in element_leaves) for clade_leaf in clade_leaves):
                    # There is a tuple wide match
                    index = [element_leaves[0] in clade_leaf for clade_leaf in clade_leaves].index(True)

                    branch_colors.append(clade_colors[index])
                else:
                    branch_colors.append("black")
        else:
            branch_colors = [tuple_color for element in tup]
    else:
        tuple_color = None
        branch_colors = ["black" for element in tup]


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
    if tuple_color:
        vcolors.append(tuple_color)
    else:
        vcolors.append("black")

    # Computing offsets
    v_offsets = [origin[1]-(v_line_length/2)+ sum(leaf_counts[:j]) + 0.5*leaf_counts[j] for j in range(len(leaf_counts))]

    for index,element in enumerate(tup): # We now cycle through again
        branch_length = float(element[1])
        horizontal_lines.append([(origin[0],v_offsets[index]),
                                 (origin[0]+branch_length,v_offsets[index])])
        hcolors.append(branch_colors[index])

        if not isinstance(element[0],str):
            new_lines = compute_segments(element[0],origin=(origin[0]+branch_length,v_offsets[index]),recursion_number=recursion_number+1,clades=clades,clade_colors=clade_colors)
            vertical_lines += new_lines[0]
            horizontal_lines += new_lines[1]
            branch_endpoints += new_lines[2]
            vcolors += new_lines[3]
            hcolors += new_lines[4]
        else:
            branch_endpoints.append((origin[0]+branch_length,v_offsets[index]))


    ### Managing unit changes
    vertical_lines = [[(u[0][0]*w_unit,u[0][1]*v_unit),(u[1][0]*w_unit,u[1][1]*v_unit)] for u in vertical_lines]
    horizontal_lines = [[(u[0][0] * w_unit, u[0][1] * v_unit), (u[1][0] * w_unit, u[1][1] * v_unit)] for u in
                      horizontal_lines]
    branch_endpoints = [(u[0]*w_unit,u[1]*v_unit) for u in branch_endpoints]

    return [vertical_lines,horizontal_lines,branch_endpoints,vcolors,hcolors]






### VISUALIZATION FUNCTIONS ###
def plot_tree(tree,
              save=True,
              labels=True,
              text_offset=0.01,
              colormode="CLADES",
              clade_point=0.0001,
              include_clade_line=True):
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
    if colormode == "CLADES":
        clds = tree.find_clades(clade_point)
        if clds:
            lines = compute_segments(tree.tree,clades=clds,clade_colors=get_random_colors(len(clds)))
        else:
            lines = compute_segments(tree.tree)
    else:
        lines = compute_segments(tree.tree)

    xs = [lines[0][i][1][0] for i in range(len(lines[0]))] + [lines[1][i][1][0] for i in range(len(lines[1]))]
    ys = [lines[0][i][1][1] for i in range(len(lines[0]))] + [lines[1][i][1][1] for i in range(len(lines[1]))]
    ax1.add_collection(LC(lines[0] + lines[1],colors=lines[3]+lines[4]))
    ## Managing spines and tick marks
    ax1.spines[:].set_visible(False)
    ax1.set_yticklabels([])
    ax1.set_yticks([])
    ax1.set_xticks([])
    ax1.set_xticklabels([])

    ## Managing Clade Line ##
    if include_clade_line:
        ax1.vlines(x=clade_point,ymin=np.amin(ys)*1.1,ymax=np.amax(ys)*1.1,color="r",ls=":")
    ## Adding arrow
    ax1.arrow(0, np.amin(ys) * (1.1), np.amax(xs) * (1.1), 0, hatch="+",
              head_width=(np.abs(np.amax(ys) - np.amin(ys)) / 50),
              head_length=np.abs(np.amax(xs) - np.amin(xs)) / 100, length_includes_head=False, fill=True,
              facecolor="k")
    ## Managing text
    ax1.set_title("Phylogenetic Tree of %s" % tree.name)
    ax1.set_xlabel("Genetic Distance Measurement")
    if labels:
        # We are actually going to add labels
        for index, label in enumerate(tree.leaves):
            ax1.text(float(lines[2][index][0]) * (1 + text_offset), lines[2][index][1], label)
    ax1.autoscale()
    plt.show()


if __name__ == '__main__':
    #log.basicConfig(level=log.DEBUG)
    #t = phylo.Tree("/home/ediggins/BioInformatics/BioPython/tree.dnd",name="test")
    t = phylo.Tree("/media/Mercury/SRR_SALMONELLA_FILES/salmonella.dnd",name="test")
    plot_tree(t,clade_point=0.00003,labels=False)