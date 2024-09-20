import matplotlib.pyplot as plt
from Bio import AlignIO
import numpy as np
from scipy.ndimage import uniform_filter1d
import seaborn as sns
import pandas as pd
import logomaker

# Read the multiple sequence alignment from a FASTA file
alignment = AlignIO.read("DAGphylum_subset_aligned_sequences.fasta", "fasta")

# Initialize lists to store gap fractions and conservation scores
gap_fractions = []
conservation_scores = []

# Function to calculate conservation score (e.g., inverse of Shannon entropy)
def calculate_conservation(column):
    # Count the frequency of each amino acid in the column, excluding gaps
    amino_acid_counts = {aa: column.count(aa) for aa in set(column) if aa != '-'}
    total_count = sum(amino_acid_counts.values())
    if total_count == 0:  # Avoid division by zero if there are only gaps
        return 0
    # Calculate the frequency of each amino acid and their respective entropy contribution
    entropy = -sum((count/total_count) * np.log(count/total_count) for count in amino_acid_counts.values())
    # Normalize the entropy (max entropy occurs with equal distribution among all amino acids)
    max_entropy = np.log(len(amino_acid_counts))
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 1
    # Conservation score is the complement of normalized entropy
    conservation_score = 1 - normalized_entropy
    return conservation_score

# Analyze each column of the alignment
for i in range(len(alignment[0])):
    column = [record.seq[i] for record in alignment]
    gap_fraction = column.count('-') / len(column)
    gap_fractions.append(gap_fraction)
    
    conservation_score = calculate_conservation(column)
    conservation_scores.append(conservation_score)

# Calculate the moving average for smoothing the line
smoothed_gaps = uniform_filter1d(gap_fractions, size=25)  # Adjust the size as needed
smoothed_conservation = uniform_filter1d(conservation_scores, size=25)  # Adjust the size as needed

# Prepare the DataFrame
data = pd.DataFrame({
    'Position': range(len(gap_fractions)),
    'Smoothed Gaps': smoothed_gaps,
    'Smoothed Conservation': smoothed_conservation
})

# Plot using Seaborn with two different y-axes
fig, ax1 = plt.subplots(figsize=(14, 6))

# Plot the 'Smoothed Gaps' on the primary y-axis (left)
sns.lineplot(data=data, x='Position', y='Smoothed Gaps', ax=ax1, color='#ff9832', label='Gap Trend', linewidth=2.5)
ax1.set_ylabel('Gap Fraction\n', color='#ff9832', fontsize=15, fontweight='bold')
ax1.tick_params(axis='y', labelcolor='#ff9832')

# Create a second y-axis for the 'Smoothed Conservation' (right)
ax2 = ax1.twinx()
sns.lineplot(data=data, x='Position', y='Smoothed Conservation', ax=ax2, color='#3299ff', label='Conservation Trend', linewidth=2.5)
ax2.set_ylabel('\nConservation Score', color='#3299ff', fontsize=15, fontweight='bold')
ax2.tick_params(axis='y', labelcolor='#3299ff')

# Remove grid
ax1.grid(False)
ax2.grid(True)

# Set the x-axis limits to the range of your positions
ax1.set_xlim(data['Position'].min(), data['Position'].max())
ax1.set_xlabel('\nPosition in MSA', color='black', fontsize=15, fontweight='bold')

plt.title('Gaps and Conservation across MSA\n', fontsize=15, fontweight='bold')
fig.tight_layout()  # Adjust the layout to fit

# Remove the legends
ax1.get_legend().remove()
ax2.get_legend().remove()

plt.show()

###

# New function for simple conservation scoring
def calculate_simple_conservation(column):
    if len(column) == column.count('-'):  # Avoid division by zero if there are only gaps
        return 0
    amino_acid_counts = {aa: column.count(aa) for aa in set(column) if aa != '-'}
    most_common_aa_count = max(amino_acid_counts.values())
    conservation_score = most_common_aa_count / len(column)
    return conservation_score

# Initialize a list for the simple conservation scores
simple_conservation_scores = []

# Calculate simple conservation scores for each column
for i in range(len(alignment[0])):
    column = [record.seq[i] for record in alignment]
    
    simple_conservation_score = calculate_simple_conservation(column)
    simple_conservation_scores.append(simple_conservation_score)

# Calculate the moving average for smoothing
smoothed_simple_conservation = uniform_filter1d(simple_conservation_scores, size=25)  # Adjust the size as needed

# Add the new scores to the DataFrame
data['Smoothed Simple Conservation'] = smoothed_simple_conservation

# Plotting for the new figure with simple conservation
fig2, ax1_2 = plt.subplots(figsize=(14, 6))

sns.lineplot(data=data, x='Position', y='Smoothed Gaps', ax=ax1_2, color='#ff9832', label='Gap Trend', linewidth=2.5)
ax1_2.set_ylabel('Gap Fraction\n', color='#ff9832', fontsize=15, fontweight='bold')
ax1_2.tick_params(axis='y', labelcolor='#ff9832')

ax2_2 = ax1_2.twinx()
sns.lineplot(data=data, x='Position', y='Smoothed Simple Conservation', ax=ax2_2, color='green', label='Conservation Trend', linewidth=2.5)
ax2_2.set_ylabel('\nConservation Score (Simple)', color='green', fontsize=15, fontweight='bold')
ax2_2.tick_params(axis='y', labelcolor='green')

ax1_2.grid(False)
ax2_2.grid(True)

ax1_2.set_xlim(data['Position'].min(), data['Position'].max())
ax1_2.set_xlabel('\nPosition in MSA', color='black', fontsize=15, fontweight='bold')

plt.title('Gaps and Conservation (Simple) across MSA\n', fontsize=15, fontweight='bold')
fig2.tight_layout()

# Remove legends
ax1_2.get_legend().remove()
ax2_2.get_legend().remove()

plt.show()

###

def create_sequence_logo(msa_file, start, end):
    # Read the MSA file
    alignment = AlignIO.read(msa_file, "fasta")

    # Extract the specified range for each sequence
    subalignments = [str(record.seq[start:end]) for record in alignment]

    # Count amino acid frequencies in each position
    position_frequencies = {i-start: {} for i in range(start, end)}
    for subalign in subalignments:
        for position, aa in enumerate(subalign):
            if aa != '-':  # Exclude gaps
                if aa not in position_frequencies[position]:
                    position_frequencies[position][aa] = 0
                position_frequencies[position][aa] += 1

    # Normalize frequencies and create DataFrame
    logo_data = pd.DataFrame(position_frequencies).fillna(0).apply(lambda x: x / len(subalignments))

    # Create the sequence logo
    logo = logomaker.Logo(logo_data.T, color_scheme='weblogo_protein')
    logo.ax.set_ylabel('Frequency')
    logo.fig.suptitle('Sequence Logo')
    
    # Remove the grid
    logo.ax.grid(False)

    # Remove x-axis tick labels
    logo.ax.set_xticks([])
    logo.ax.set_xticklabels([])
    
    # Show the plot
    plt.show()

# C-term Extension All Sequences
create_sequence_logo('DAG_all_aligned_sequences.fasta', 2070, 2160)

# Auto-Proteolysis Site All Sequences
create_sequence_logo('DAG_all_aligned_sequences.fasta', 1784, 1800)

# C-term Extension Phylum
create_sequence_logo('DAGphylum_subset_aligned_sequences.fasta', 795, 820)

# Auto-Proteolysis Site Phylum
create_sequence_logo('DAGphylum_subset_aligned_sequences.fasta', 737, 755)
