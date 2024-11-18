import os
import json
import typing

TEMPFILE = 'temp_pdb'
output_file = './output.txt'
directory_path = '/research/jagodzinski/DATA/mutants/1hhp/1/2'

class Atom:
    def __init__(self, name : str, x : float, y : float, z : float):
        self.name = name
        self.x = x
        self.y= y
        self.z = z

    def midpoint(self, other : Type[Atom]):
        # TODO
        ...

    def distance(self, other : Atom):
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
                        content = json.load(input_file)['pdb_data']['pdb']
                        temp_pdb.write(content)
                        temp_pdb.seek(0, 0)
                        os.system(f"./hbplus {TEMPFILE}.pdb")

                        output_file.write(parse_hb_file(f"{TEMPFILE}.hb2"))
                    
                    output_file.write("\n" + "-" * 50 + "\n")
            
        print(f"File created and updated successfully at: {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}\n")


def parse_hb_file(file):
    with open(file, 'r') as hb_data:

        return hb_data.read()





if __name__ == "__main__":
    main()