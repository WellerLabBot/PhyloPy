"""

    Python Package for Managing SRA Toolkit actions with python.
            Written By: Eliza Diggins

"""
import numpy as np
import pathlib as pt
import pandas as pd
import logging as log
import os
from tqdm import tqdm

### BASE FUNCTIONS ###
def gzip(directory,filetype=".fastq"):
    """
    G zips all files of type filetype in the directory.
    :param directory: The directory to zip into
    :param filetype: The filetype to zip
    :return: None
    """
    print("BioPython:Tools:BioPy_SRA:gzip: Attempting to zip all %s files in %s."%(filetype,directory))

    try:
        os.chdir(directory)
    except Exception:
        print("BioPython:Tools:BioPy_SRA:gzip:ERROR failed to find directory %s. Please try again."%directory)

    files = [file for file in os.listdir() if file[-len(filetype):] == filetype]

    print("BioPython:Tools:BioPy_SRA:gzip: Found %s files with extension %s."%(len(files),filetype))
    fails = []
    for file in tqdm(files,desc="Zipping files."):
        tqdm.write("BioPython:Tools:BioPy_SRA:gzip: Zipping %s."%file)
        try:
            os.system("gzip %s 1>/dev/null 2>/dev/null"%file)
        except Exception:
            tqdm.write("BioPython:Tools:BioPy_SRA:gzip:WARNING: zip failed on %s."%file)
            fails.append(file)

    print("BioPython:Tools:BioPy_SRA:gzip: Finished zipping %s. Successfully zipped %s of %s (%% %s)."%(directory,len(files)-len(fails),len(files),100*(1-(len(fails)/len(files)))))
    return True


### FUNCTIONS ###
def download_fasterq(inp,output_location,options=None,echo=True):
    """
    Downloads the SRA files and compiles into fasterq.
    :param inp: The inp, either an ascension or a list of ascensions, or a file of ascensions.
    :param output_location: The location to store the output files
    :return: True if passed, False if not.
    """
    log.debug("BioPython:Tools:BioPy_SRA:download_fasterq:DEBUG: Attempting to download data from %s to %s."%(inp,output_location))

    ### Sanitizing inp ###
    # Managing inp
    inp_TYPE = None # This will mark the inp type when we know it.
    if isinstance(inp,(str)): # The inp is a string. It could be either a filepath or an ascension.
        # Checking if file.
        if os.path.isfile(inp):
            log.debug("BioPython:Tools:BioPy_SRA:download_fasterq:DEBUG: %s is identified as a file."%inp)

            # Reading the file.
            log.debug("BioPython:Tools:BioPy_SRA:download_fasterq:DEBUG: Reading %s"%inp)

            try:
                with open(inp,"r+") as file:
                    ascs = file.read().splitlines()
            except Exception:
                log.error("BioPython:Tools:BioPy_SRA:download_fasterq:ERROR: Failed to read input file %s."%inp)
                return False

        elif inp[:3] == "SRR":
            log.debug("BioPython:Tools:BioPy_SRA:download_fasterq:DEBUG: %s is identified as an ascension."%inp)
            ascs = [inp]
        else:
            log.error("BioPython:Tools:BioPy_SRA:download_fasterq:ERROR: %s was of type %s, but was neither a file nor an ascension. Returning False."%(inp,type(inp)))
            return False
    elif isinstance(inp,(list)):
        # The inp was a list, we need to check that they are all valid.
        if any(asc[:3] != "SRR" for asc in inp):
            # The inps did not all have the right prefix/
            log.error(
                "BioPython:Tools:BioPy_SRA:download_fasterq:ERROR: Ascension list %s had invalid ascensions. Returning False." % (inp))
            return False
        else:
            log.debug("BioPython:Tools:BioPy_SRA:download_fasterq:DEBUG: %s is identified as an list."%inp)
            ascs = inp
    else:
        log.error(
            "BioPython:Tools:BioPy_SRA:download_fasterq:ERROR: Input %s was not a valid input type. See documentation." % (
                inp))
        return False

    # Managing output_location
    if isinstance(output_location,str): # the type is correct
        if not os.path.isdir(output_location):
            # The path doesn't already exist
            log.warning('BioPython:Tools:BioPy_SRA:download_fasterq:WARNING: Failed to find %s. Attempting to create the directory.'%output_location)

            try:
                os.mkdir(output_location)
            except:
                log.error("BioPython:Tools:BioPy_SRA:download_fasterq:ERROR: output_location: %s could not be made. Returning False."%output_location)
                return False
    else:
        log.error(
            "BioPython:Tools:BioPy_SRA:download_fasterq:ERROR: output_location: %s was not a valid input. Returning False." % output_location)
        return False

    # Managing Options
    if not options: # The options are being set to default. #TODO: git and all else.
        options = {
            "kwargs":["--concatenate-reads","-t '/media/Mercury/SRR_TMP_CACHE'"]
        }

    # building kwarg string
    kwgs = ""
    for kwg in options["kwargs"]:
        kwgs += kwg
        kwgs += " "


    ### FULL RUN ###

    log.info("BioPython:Tools:BioPy_SRA:download_fasterq:INFO: Attempting to run FATERQ on %s"%ascs)
    if echo:
        print("BioPython:Tools:BioPy_SRA:download_fasterq:INFO: Attempting to run FATERQ. Ascensions:\n\n")
        print("##########################################  ASCENSIONS  ##############################################")
        for asc in ascs:
            print(asc)

        print("######################################################################################################")

    for asc in tqdm(ascs,desc="Downloading Ascensions"):
        # Attempting downloads.
        log.debug("BioPython:Tools:BioPy_SRA:download_fasterq:DEBUG: fasterq-dump %s -O '%s' %s"%(asc,output_location,kwgs))
        tqdm.write("Downloading %s...   "%asc)
        os.system("fasterq-dump %s -O '%s' %s 2>/dev/null"%(asc,output_location,kwgs))

    if echo:
        print("BioPython:Tools:BioPy_SRA:download_fasterq:INFO: Completed downloads...")

def generate_mashtree(directory,output_name,echo=True):
    """
    Generates a mashtree phylogenetic tree from the files at the given directory.
    :param directory: The directory to use.
    :return: None
    """
    # Intro Logging
    log.info("BioPython:Tools:BioPy_SRA:generate_mashtree:INFO: Attempting to run FATERQ on %s" % directory)

    os.chdir(directory)
    # Printing to console.
    if echo:
        print("BioPython:Tools:BioPy_SRA:generate_mashtree:INFO: Attempting to build Mashtree File from %s:\n\n"%directory)
        print("BioPython:Tools:BioPy_SRA:generate_mashtree:INFO: Locating Files...")
        print("BioPython:Tools:BioPy_SRA:generate_mashtree:INFO: Found %s files."%len([i for i in os.listdir() if i[-6:]==".fastq"]))
        print("####################################  .fastq  FILES  ##############################################")
        for file in os.listdir():
            if file[-6:]==".fastq":
                print(file)

        print("######################################################################################################")
        print("BioPython:Tools:BioPy_SRA:generate_mashtree:INFO: Generating mashtree:")

        os.system("export PATH=$HOME/bin:$PATH")
        os.system("export PERL5LIB=$PERL5LIB:$HOME/lib/perl5")
        os.system("mashtree *.fastq > %s"%output_name)

        print("BioPython:Tools:BioPy_SRA:generate_mashtree:INFO: Saved %s in %s."%(output_name,directory))
if __name__ == '__main__':
    gzip("/media/Mercury/SRR_SALMONELLA_FILES")