from io import TextIOWrapper
import os
import json
import math

from typing import Self, Dict , Tuple

TEMPFILE = 'temp_pdb'
output_file = './output.txt'
# directory_path = '/research/jagodzinski/DATA/mutants/1hhp/1/2'
directory_path = 'data'
PDB_Dict = Dict[Tuple[str, str], Tuple[float, float, float]]

def main():
    create_and_write_file(output_file=output_file, directory_path=directory_path)

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
def create_and_write_file(output_file : str, directory_path : str):

    try:
        # Create and open output file
        with open(output_file, 'w') as output_file:

            # Make a temporary file to put pdb data from json
            with open(f"./{TEMPFILE}.pdb", 'w') as temp_pdb:

                # List all files in the specified directory
                for item in os.listdir(directory_path):
                    ins_loc1, ins_typ1, ins_loc2, ins_typ2 = parse_mutation_location(item)
                    item_path = os.path.join(directory_path, item)
                    with open(item_path, 'r') as input_file:
                        # get pdb data from json file and put it into our temp file
                        content = json.load(input_file)['pdb_data']['pdb']
                        # content = input_file.read()
                        temp_pdb.write(content)
                        temp_pdb.seek(0, 0)
                        
                        pdb_dict = build_dict(content)
                        # calculate hydrogen locations
                        os.system(f"./hbplus.exe {TEMPFILE}.pdb")

                        # get atoms that represent mutation points
                        atom1 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc1)
                        atom2 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc2)

                        # get distances to mutation points
                        distances = parse_hb_file(f"{TEMPFILE}.hb2", mutation1=atom1, mutation2=atom2, pdb_dict=pdb_dict)
                        line = " ".join([ins_loc1, ins_typ1, ins_loc2, ins_typ2, distances])
                        output_file.write(line)
                    
    except Exception as e:
        print(f"An error occurred: {e}\n")
        exit()

def find_atom(pdb_dict : PDB_Dict , atom_name : str, residue_num : str) -> Atom:
    pos = pdb_dict[(residue_num, atom_name)]
    return Atom(atom_name, *pos)
    

def build_dict(pdb_lines) -> PDB_Dict:

    dict = {}
    for line in pdb_lines.splitlines():
        if line.startswith("ATOM"):
            residue_number = line[22:26].strip()
            atom_name = line[12:16].strip()
            record = (
                float(line[30:38].strip()),
                float(line[38:46].strip()),
                float(line[46:54].strip()),
            )
            
            key = (residue_number, atom_name)
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
        data = hblus_data.readlines()[9:]
        distances = []
        for line in data:
            donor_atom_num, donor_atom_name = line[1:5].lstrip('0'), line[9:13].strip()
            acceptor_atom_num, acceptor_atom_name = line[15:19].lstrip('0'), line[23:27].strip()

            donor = find_atom(pdb_dict=pdb_dict, atom_name=donor_atom_name, residue_num=donor_atom_num)
            acceptor = find_atom(pdb_dict=pdb_dict, atom_name=acceptor_atom_name, residue_num=acceptor_atom_num)
            
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
    return mutation[9:16].split(sep="_")


if __name__ == "__main__":
    main()