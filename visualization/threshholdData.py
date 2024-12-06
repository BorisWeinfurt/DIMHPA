import csv
import traceback

residue_length = 101
residue_numbers = [str(i) for i in range(1, residue_length + 1)]
amino_acids = [
    'A',  # Alanine
    'R',  # Arginine
    'N',  # Asparagine
    'D',  # Aspartic acid
    'C',  # Cysteine
    'E',  # Glutamic acid
    'Q',  # Glutamine
    'G',  # Glycine
    'H',  # Histidine
    'I',  # Isoleucine
    'L',  # Leucine
    'K',  # Lysine
    'M',  # Methionine
    'F',  # Phenylalanine
    'P',  # Proline
    'S',  # Serine
    'T',  # Threonine
    'W',  # Tryptophan
    'Y',  # Tyrosine
    'V'   # Valine
]

def process_data_for_heatmap(input_file, output_name, threshold, is_residue_num):
    """
    Reads lines from input_file, processes each line, and writes the result to output_file in CSV format.
    """
    my_lst = residue_numbers if is_residue_num else amino_acids
    above_threshold = {key: {k: 0 for k in my_lst} for key in my_lst}
    below_threshold = {key: {k: 0 for k in my_lst} for key in my_lst}
    try:
        with open(input_file, 'r') as infile:
            lines = infile.readlines()
            for index, line in enumerate(lines):
                split_line = line.split(" ")
                split_line[-1] = split_line[-1].strip("\n")
                
                loc1, typ1, loc2, typ2 = split_line[0:4]
                for distance in split_line[4:]:
                    try:
                        # Attempt to convert the distance to a float
                        float_distance = float(distance)
                        if is_residue_num:
                            if float_distance >= threshold:
                                above_threshold[loc1][loc2] += 1
                            else:
                                below_threshold[loc1][loc2] += 1
                        else:
                            if float_distance >= threshold:
                                above_threshold[typ1][typ2] += 1
                            else:
                                below_threshold[typ1][typ2] += 1
                    except ValueError:
                        print(f"Error converting {index} distance: \n'{distance}'\n to float in line:\n {line}\n Split line:\n {split_line}")
                        exit()
                            
        write_square_dict_to_csv(above_threshold, output_name + "_above_thresh.csv", my_lst)
        write_square_dict_to_csv(below_threshold, output_name + "_below_thresh.csv", my_lst)
    except Exception as e:
        print("EERRRROORR:", e)
    
def write_square_dict_to_csv(matrix, filename, keys):
    """
    Writes a dictionary of dictionaries to a CSV file.
    The outer dictionary keys form the rows, and the inner dictionary keys form the columns.
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow(['Second Insertion'] + keys)
        
        # Write the rows
        for row_key in keys:
            row = [row_key] + [matrix[row_key][col_key] for col_key in keys]
            writer.writerow(row)     
        

def process_data_for_histogram(input_file, output_name):
    """
    Reads lines from input_file, processes each line, and writes the result to output_file in CSV format.
    """
    ordered_acids = ['A', 'G', 'S', 'N', 'D', 'C', 'P', 'T', 'Q', 'E', 'H', 'V', 'R', 'I', 'L', 'K', 'M', 'F', 'W', 'Y']
    num_buckets = 35
    buckets = range(num_buckets)
    my_dict = {"".join(sorted([i, j])): {key : 0 for key in buckets} for i in ordered_acids for j in ordered_acids}


    try:
        with open(input_file, 'r') as infile:
            for index, line in enumerate(infile):
                split_line = line.split(" ")
                split_line[-1] = split_line[-1].strip("\n")
                
                loc1, typ1, loc2, typ2 = split_line[0:4]
                for distance in split_line[4:]:
                    try:
                        # Attempt to convert the distance to a float
                        float_distance = float(distance)
                        bucket = min(num_buckets - 1, int(float_distance))
                        key = "".join(sorted([typ1, typ2]))
                        
                        my_dict[key][bucket] += 1
                    except ValueError:
                        print(f"Error converting {index} distance: \n'{distance}'\n to float in line:\n {line}\n Split line:\n {split_line}")
        print("finished")              
        write_histo_dict_to_csv(my_dict, output_name + ".csv", buckets)
    except Exception as e:
        traceback.print_exc()
        print("EERRRROORR:", e, "\n", line)
        
def write_histo_dict_to_csv(histo_dict, output_file, bucket_list):
    """
    Writes a histogram dictionary to a CSV file.
    
    Args:
        histo_dict (dict): The histogram dictionary to write. 
                           Outer keys map to inner dictionaries of bucket counts.
        output_file (str): The file path to save the CSV.
        bucket_list (list or range): A list or range of bucket labels for the header.
    """
    try:
        with open(output_file, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the header row
            header = ["Key"] + [f"{bucket}" for bucket in bucket_list]
            writer.writerow(header)

            # Write each row of data
            for key, bucket_counts in histo_dict.items():
                row = [key] + [bucket_counts[bucket] for bucket in bucket_list]
                writer.writerow(row)

    except Exception as e:
        print(f"Error writing histogram data to {output_file}: {e}")

if __name__ == "__main__":
    input_path = "full_output.txt"
    output_path = "1hhp"
    threshold = 15

    # process_data_for_heatmap(input_path, output_path + "_residue_type", threshold=threshold, is_residue_num=False)
    # process_data_for_heatmap(input_path, output_path + "_residue_location", threshold=threshold, is_residue_num=True)
    process_data_for_histogram(input_path, output_path + "_histogram")