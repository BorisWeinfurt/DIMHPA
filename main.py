from io import TextIOWrapper
import os
import json
from typing import Self

TEMPFILE = 'temp_pdb'
output_file = './output.txt'
directory_path = '/research/jagodzinski/DATA/mutants/1hhp/1/2'

class Atom:
    def __init__(self, name : str, x : float, y : float, z : float):
        self.name = name
        self.x = x
        self.y= y
        self.z = z

    def midpoint(self, other : Self):
        # TODO
        ...

    def distance(self, other : Self):
        # TODO
        ...


def main():
    create_and_write_file(output_file=output_file, directory_path=directory_path)

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
                    item_path = os.path.join(directory_path, item)
                    with open(item_path, 'r') as input_file:

                        # get pdb data from json file
                        content = json.load(input_file)['pdb_data']['pdb']
                        temp_pdb.write(content)
                        temp_pdb.seek(0, 0)

                        # calculate hydrogen locations
                        os.system(f"./hbplus {TEMPFILE}.pdb > /dev/null")

                        # get atoms that represent mutation points
                        temp_pdb.seek(0, 0)
                        atom1 = None
                        atom2 = None

                        # get distances to mutation points
                        distances = parse_hb_file(f"{TEMPFILE}.hb2",atom1,atom2)
                        output_file.write(distances)
                    
                    output_file.write("\n" + "-" * 50 + "\n")
            
    except Exception as e:
        print(f"An error occurred: {e}\n")
        exit()

"""
Find a specific atom in a pdb file based on an atom name and residue num

:param file: open file buffer to the target pdb file
:param atom_name: The atom type 4 letter identifier to search for
:param resigue_num: The residue number to search for
"""
def find_atom(file : TextIOWrapper, atom_name : str, residue_num : int) -> Atom:
    ...


"""
Look through each hydrogen bond outputted by HBPLUS and calculate the minimum 
distance to a mutation point. 

:param file: open file buffer to the target HBPLUS output file
:param mutation1: Atom representing the first mutation
:param mutation2: Atom representing the second mutation
"""
def parse_hb_file(file : TextIOWrapper, mutation1 : Atom, mutation2 : Atom):
    with open(file, 'r') as hb_data:

        return hb_data.read()





if __name__ == "__main__":
    main()