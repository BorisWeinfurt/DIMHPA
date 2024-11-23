from io import TextIOWrapper
import os
import json
import re
import math

from typing import Self, Dict , Tuple

TEMPFILE = 'temp_pdb'
output_file = './output.txt'
directory_path = '/research/jagodzinski/DATA/mutants/1hhp/100'
# directory_path = 'data'
PDB_Dict = Dict[Tuple[str, str], Tuple[float, float, float]]

def main():
    create_and_write_file(output_file=output_file, directory_path=directory_path, tempfile=TEMPFILE)

class Atom:
    def __init__(self, name : str, x : float, y : float, z : float):
        self.name = name
        self.x = x
        self.y= y
        self.z = z
    
    def __repr__(self):
        return f"Atom(name={self.name}, x={self.x}, y={self.y}, z={self.z})"

    def midpoint(self, other: Self):
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
Create and open a file at a specific path, list every file in a specified directory,
and append content to the file.

:param file_path: The full path to the file to create and write.
:param directory_path: The directory to look for files.
:param append_content: The content to append to the file.
"""
def create_and_write_file(output_file : str, directory_path : str, tempfile=str):

    try:
        # Create and open output file
        open(output_file, 'w').close()
        with open(output_file, 'a') as output_file:

            # Make a temporary file to put pdb data from json
            with open(f"./{tempfile}.pdb", 'w') as temp_pdb:

                # List all files in the specified directory
                for (dirpath, dirnames, filenames) in os.walk(directory_path):
                    print(f"analyzing: dirpath-{dirpath}")
                    for file in filenames:
                        item = dirpath + "/" + file

                        ins_loc1, ins_typ1, ins_loc2, ins_typ2 = parse_mutation_location(item)

                        item_path = os.path.join(directory_path, item)
                        with open(item_path, 'r') as input_file:


                            # get pdb data from json file and put it into our temp file
                            content = json.load(input_file)['pdb_data']['pdb']
                            # content = input_file.read()
                            temp_pdb.write(content)
                            temp_pdb.seek(0, 0)
                            

                            # calculate hydrogen locations
                            os.system(f"./hbplus {tempfile}.pdb -o > err")
                            pdb_dict = build_dict("./" + tempfile + ".h")
                            # get atoms that represent mutation points
                            atom1 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc1, residue_name=one_to_three(ins_typ1))
                            atom2 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc2, residue_name=one_to_three(ins_typ2))

                            # get distances to mutation points
                            distances = parse_hb_file(f"{tempfile}.hb2", mutation1=atom1, mutation2=atom2, pdb_dict=pdb_dict)
                            line = " ".join([ins_loc1, ins_typ1, ins_loc2, ins_typ2, distances]) + "\n"
                            output_file.write(line)
                    
    except Exception as e:
        print(f"An error occurred: {e}\n")
        exit()

def find_atom(pdb_dict : PDB_Dict , atom_name : str, residue_num : str, residue_name : str) -> Atom:
    pos = pdb_dict[(residue_num, atom_name, residue_name)]
    return Atom(atom_name, *pos)
    

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
            split = re.split('\s+', line)
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


def one_to_three(one_letter_code):
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
    return mapping[one_letter_code.upper()]

if __name__ == "__main__":
    main()