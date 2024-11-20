from io import TextIOWrapper
import os
import json
import math

from typing import Self, Dict , Tuple

TEMPFILE = 'temp_pdb'
output_file = './output.txt'
directory_path = '/research/jagodzinski/DATA/mutants/1hhp/1/2'
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

    def midpoint(self, other: 'Atom'):
        """
        Calculate the midpoint between this atom and another atom.
        Returns a new Atom representing the midpoint.
        """
        mid_x = (self.x + other.x) / 2
        mid_y = (self.y + other.y) / 2
        mid_z = (self.z + other.z) / 2
        return (mid_x, mid_y, mid_z)

    def distance(self, other : Self):
        # TODO
        ...


def create_and_write_file(output_file : str, directory_path : str):
    """
    Create and open a file at a specific path, list every file in a specified directory,
    and append content to the file.
    
    :param file_path: The full path to the file to create and write.
    :param directory_path: The directory to look for files.
    :param append_content: The content to append to the file.
    """
    try:
        # Create and open output file
        with open(output_file, 'a') as output_file:

            # Make a temporary file to put pdb data from json
            with open(f"./{TEMPFILE}.pdb", 'w') as temp_pdb:

                # List all files in the specified directory
                for item in os.listdir(directory_path):
                    ins_loc1, ins_typ1, ins_loc2, ins_typ2 = parse_mutation_location(item)
                    item_path = os.path.join(directory_path, item)

                    with open(item_path, 'r') as input_file:

                        # get pdb data from json file and put it into our temp file
                        content = json.load(input_file)['pdb_data']['pdb']
                        temp_pdb.write(content)
                        temp_pdb.seek(0, 0)

                        pdb_dict = build_dict(content)

                        # calculate hydrogen locations
                        os.system(f"./hbplus {TEMPFILE}.pdb > /dev/null")

                        # get atoms that represent mutation points
                        atom1 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc1)
                        atom2 = find_atom(pdb_dict=pdb_dict, atom_name='CA', residue_num=ins_loc2)

                        print(atom1, atom2)
                        # get distances to mutation points
                        distances = parse_hb_file(f"{TEMPFILE}.hb2", atom1, atom2)
                        output_file.write(distances)
                    
                    output_file.write("\n" + "-" * 50 + "\n")
            
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

:param file: open file buffer to the target HBPLUS output file
:param mutation1: Atom representing the first mutation
:param mutation2: Atom representing the second mutation
"""
def parse_hb_file(file : TextIOWrapper, mutation1 : Atom, mutation2 : Atom):
    # TODO
    with open(file, 'r') as hb_data:

        return hb_data.read()

"""
Given a JSON file from the research directory determing the locations and mutations being inserted

:param mutation: name of the file currently being processed"""
def parse_mutation_location(mutation : str) -> tuple[str, str, str, str]:
    return mutation[9:16].split(sep="_")


if __name__ == "__main__":
    main()