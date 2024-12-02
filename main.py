import glob
from io import TextIOWrapper
from multiprocessing import Process
import os
import json
import re
import math

import shutil
import sys
from typing import Dict , Tuple

PDB_Dict = Dict[Tuple[str, str], Tuple[float, float, float]]

"""
Map the one letter abbreviations of amino acids to 3 letter ones
"""
mapping = {'A':'ALA',
               'R':'ARG',
               'N':'ASN',
               'D':'ASP',
               'C':'CYS',
               'Q':'GLN',
               'E':'GLU',
               'G':'GLY',
               'H':'HIS',
               'I':'ILE',
               'L':'LEU',
               'K':'LYS',
               'M':'MET',
               'F':'PHE',
               'P':'PRO',
               'S':'SER',
               'T':'THR',
               'W':'TRP',
               'Y':'TYR',
               'V':'VAL',}


class Atom:
    def __init__(self, name : str, x : float, y : float, z : float):
        self.name = name
        self.x = x
        self.y= y
        self.z = z
    
    def __repr__(self):
        return f"Atom(name={self.name}, x={self.x}, y={self.y}, z={self.z})"

    def midpoint(self, other):
        """
        Calculate the midpoint between this atom and another atom.
        Returns a new Atom representing the midpoint.
        """
        mid_x = (self.x + other.x) / 2
        mid_y = (self.y + other.y) / 2
        mid_z = (self.z + other.z) / 2
        return (mid_x, mid_y, mid_z)

    def distance(self, other : Tuple[float,float,float]) -> float:
        return math.sqrt(((self.x - other[0])**2 + (self.y - other[1])**2 + (self.z - other[2])**2))
        

"""
Analyzes a given directory from the /research/jagodzinski/DATA/mutants/*/* path
this means that it should be a directory of directories that contain pdbs stored as json files

Each pdb is obtained from json, run through hbplus, and then distances are analyzed
via the hbplus output files

:param output_file: the name of the output file
:param tempfile: the name of the temporary file
:param paths_to_analyze: a list of paths of filenames this process should analyze
"""
def anaylze_file(output_file : str, tempfile : str, paths_to_analyze : list):
    
    try:
        # Create and open output file
        with open(output_file, 'a') as output_file:

            # Make a temporary file to put pdb data from json
            with open(f"{tempfile}.pdb", 'w') as temp_pdb:

                # List all files in the specified directory
                for file_path in paths_to_analyze:
                    ins_loc1, ins_typ1, ins_loc2, ins_typ2 = parse_mutation_location(file_path)

                    with open("../" + file_path, 'r') as input_file:

                        # get pdb data from json file and put it into our temp file
                        content = json.load(input_file)['pdb_data']['pdb']
                        # content = input_file.read()
                        temp_pdb.write(content)
                        temp_pdb.seek(0, 0)
                        
                        # calculate hydrogen locations
                        os.system(f"../hbplus {tempfile}.pdb -o > err")

                        # use .h file instead of original pdb to both
                        # 1. account for possible uncertainty/duplicates in original pdb
                        # 2. account for hydrogens that are not in the original pdb
                        pdb_dict = build_dict(tempfile + ".h")
                        # get atoms that represent mutation points
                        atom1 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc1, residue_name=mapping[ins_typ1])
                        atom2 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc2, residue_name=mapping[ins_typ2])

                        # get distances to mutation points
                        distances = parse_hb_file(f"{tempfile}.hb2", mutation1=atom1, mutation2=atom2, pdb_dict=pdb_dict)
                        line = " ".join([ins_loc1, ins_typ1, ins_loc2, ins_typ2, distances]) + "\n"
                        output_file.write(line)
    except Exception as e:
        print(f"An error occurred: {e}\n")

"""
Find an atom from the pdb dictionary given the parameters
"""
def find_atom(pdb_dict : PDB_Dict , atom_name : str, residue_num : str, residue_name : str) -> Atom:
    pos = pdb_dict[(residue_num, atom_name, residue_name)]
    return Atom(atom_name, *pos)
    
"""
Build a dictionary of atoms from a pdb file
NOTE: this is inteded to be used with the hblus output file NOT regular pdbs
"""
def build_dict(pdb_file) -> PDB_Dict:

    dict = {}
    with open(pdb_file, 'r') as pdb:
        for line in pdb.readlines():
            if line.startswith("ATOM"):
                split = re.split('\s+', line)
                residue_number = split[5].strip()
                atom_name = split[2].strip()
                residue_name = split[3].strip()

                record = (
                    float(split[6].strip()),
                    float(split[7].strip()),
                    float(split[8].strip()),
                )
                key = (residue_number, atom_name, residue_name)
                dict[key] = record
    return dict

"""
Look through each hydrogen bond outputted by HBPLUS and calculate the minimum 
distance to a mutation point. 

:param file: file name of hbplus result
:param mutation1: Atom representing the first mutation
:param mutation2: Atom representing the second mutation
:param pdb_dict: Dictionary of the current pdb file
"""
def parse_hb_file(hb_file : str, mutation1 : Atom, mutation2 : Atom, pdb_dict : PDB_Dict):
    with open(hb_file, 'r') as hblus_data:
        # skip the headers
        data = hblus_data.readlines()[8:]
        distances = []
        for line in data:
            split = re.split('\\s+', line)
            donor_atom_res_num, donor_atom_name, donor_res_typ = split[0][1:5].lstrip('0'), split[1], split[0][6:10]
            acceptor_atom_res_num, acceptor_atom_name, acceptor_res_typ = split[2][1:5].lstrip('0'), split[3], split[2][6:10]

            donor = find_atom(pdb_dict=pdb_dict, atom_name=donor_atom_name, residue_num=donor_atom_res_num, residue_name=donor_res_typ)
            acceptor = find_atom(pdb_dict=pdb_dict, atom_name=acceptor_atom_name, residue_num=acceptor_atom_res_num, residue_name=acceptor_res_typ)
            
            hydrogen_position = donor.midpoint(acceptor)
            
            min_distance = min(
                mutation1.distance(hydrogen_position),
                mutation2.distance(hydrogen_position),
                )
            distances.append(min_distance)
            
        return " ".join(map(str, distances))
            

"""
Given a JSON file from the research directory determing the locations and mutations being inserted

:param mutation: name of the file currently being processed"""
def parse_mutation_location(mutation : str) -> tuple[str, str, str, str]:
    return re.split("[_.]", mutation)[2:6]


"""
Gives each process a file to do read/write to and splits up the input files based on how many directories there are
params: 
id - the id of the process
file_list - the paths to the files this process should analyze
num_procs - the number of processes working on the task
"""

def wrapper(id,file_list,result_directory):
    # generate file names for work
    temp_file = "temp_file_" + str(id)
    out_file = 'outfile_' + str(id)
    
    #this is our final list with what we should need
    os.chdir(result_directory)
    anaylze_file(output_file=out_file, tempfile=temp_file,paths_to_analyze=file_list)
    
    os.remove(temp_file + ".h")
    os.remove(temp_file + ".hb2")
    os.remove(temp_file + ".pdb")


def get_file_paths(directory):
    """Recursively get all file paths in the inner directories."""
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            file_paths.append(os.path.join(root, file))
    return file_paths

def split_files(file_paths, num_processes):
    """Split the file paths into chunks for each process."""
    chunk_size = len(file_paths) // num_processes
    remainder = len(file_paths) % num_processes
    chunks = []
    start = 0
    
    for i in range(num_processes):
        end = start + chunk_size + (1 if i < remainder else 0)
        chunks.append(file_paths[start:end])
        start = end
    return chunks

def combine_data(data_path, output_file):

    with open(output_file, 'wb') as outfile:
        for filename in glob.glob(data_path + ' /*'):

            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

if __name__ == '__main__':
    # Just boot up however many processes we want and have them start the function
    # the math in the wrapper covers generating names for all the temp and output files as well as how many folders each process needs to examine
    # all the values we should have to touch are the number of processes and the protein protein_length
    # as well as any paths used anywhere
    if len(sys.argv) < 2:
        print("Usage: python script.py <directory>")
        exit(1)
    
    num_procs = 7
    p_list = []
    
    mutant = sys.argv[1]
    directory_path = '/research/jagodzinski/DATA/mutants/'
    # directory_path = ''
    
    all_file_paths = get_file_paths(directory=directory_path+mutant)
    chunks = split_files(file_paths=all_file_paths, num_processes=num_procs)
    
    result_directory = mutant + "_results"
    os.makedirs(result_directory, exist_ok=True)
    
    # run each process on the a chunk of the files
    for i in range (0,num_procs):
        p = Process(target=wrapper, args=(i,chunks[i], result_directory))
        p.start()
        p_list.append(p)
        
    # wait for all processes to finish
    for p in p_list:
        p.join()
        
    # remove junk files
    os.remove(result_directory + "/err")
    os.remove(result_directory + "/hbdebug.dat")
    
    # combine data into a single file
    combine_data(data_path=result_directory, output_file = result_directory + '.txt')
    
    shutil.rmtree(result_directory)
    
    