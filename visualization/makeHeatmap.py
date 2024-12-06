import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# List of CSV file paths
csv_files = ["1hhp_residue_location_above_thresh.csv", "1hhp_residue_location_below_thresh.csv", "1hhp_residue_type_above_thresh.csv", "1hhp_residue_type_below_thresh.csv"]  # Replace with your file names

# Calculate the grid dimensions for square plots
n_files = len(csv_files)
rows = int(n_files**0.5)
cols = rows if rows * rows >= n_files else rows + 1

# Create a figure with subplots
fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))  # Square aspect for each heatmap
axes = axes.flatten()  # Flatten to iterate easily over the axes
# Loop through CSV files and plot each as a heatmap
for i, csv_file in enumerate(csv_files):
    matrix = pd.read_csv(csv_file, index_col=0)  # Load CSV with row labels
    sns.heatmap(matrix, annot=False, fmt=".2f", cmap="viridis", ax=axes[i])
    axes[i].set_title(csv_file)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()
