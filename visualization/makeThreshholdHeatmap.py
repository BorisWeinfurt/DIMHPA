import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def generate_heatmap(csv_file, output_image_name):
    """
    Generate a heatmap from a CSV file, transposing the matrix and then flipping the y-axis.
    The x-axis corresponds to the columns, and the y-axis is the rows flipped after transposition.

    Args:
    - csv_file (str): Path to the CSV file.
    - output_image_name (str): Path to save the generated heatmap image.
    """
    try:
        # Read the CSV file
        df = pd.read_csv(csv_file, index_col=0)
        
        # Transpose the matrix (swap rows and columns)
        df_transposed = df.T
        
        # Flip the y-axis of the transposed matrix (reverse the rows)
        df_flipped_y = df_transposed.iloc[::-1]

        # Create the heatmap
        plt.figure(figsize=(20, 6))
        sns.heatmap(df_flipped_y, annot=False, cmap="viridis", cbar=True)

        # Label axes and title
        plt.xlabel("Residue Pair")
        plt.ylabel("Distance to mutation")
        plt.title("Residue Pairs Distance Histogram")
        
        # Save the heatmap
        plt.savefig(output_image_name)
        plt.close()
        
        print(f"Heatmap saved as {output_image_name}")
    except Exception as e:
        print("Error generating heatmap:", e)



if __name__ == "__main__":
    generate_heatmap("1hhp_histogram.csv", "test.png")