import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import edlib

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list


color_map = {
    'heavy': '#FA003F',
    'light': '#00916E'
}

# Define VDJ colors
color_genes = {
    "V": "#53B3CB", 
    "D": "#F9C22E",
    "J": "#C179B9" 
}

# The assigned clusters
cluster_map = {
  '08D18': 1,
  '08O09': 1,
  '03F18': 1,
  '08H10': 1,
  '08N23': 1,
  '08K19': 1,
  '08M13': 1,
  '03L02': 1,
  '08J10': 1,
  '10H18': 1,
  '08O03': 1,
  '05K07': 1,
  '09I10': 1,
  '05N02': 2,
  '05B17': 2,
  '08F04': 2,
  '05D08': 2,
  '05D14': 2,
  '05C11': 2,
  '05M13': 2
}

# Get script path
file_path  = os.path.split(__file__)[0]
file_path = os.path.abspath(os.path.join(file_path, '..'))

# Create the output path
output_path = os.path.join(file_path, 'results')
if not os.path.exists(output_path):
  os.mkdir(output_path)
  
# Load the dataset
df = pd.read_csv(os.path.join(file_path, 'data', 'antibodies_igblast.tsv'), sep='\t')

# Define the BLOSUM62 matrix
BLOSUM62 = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1,
    ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2,
    ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
    ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0,
    ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3,
    ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
    ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0,
    ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3,
    ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
    ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1,
    ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1,
    ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
    ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1,
    ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1,
    ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
    ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2,
    ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1,
    ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
    ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1,
    ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3,
    ('E', 'Y'): -2, ('E', 'V'): -2,
    ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3,
    ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3,
    ('G', 'V'): -3,
    ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1,
    ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
    ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3,
    ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
    ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2,
    ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
    ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1,
    ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
    ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1,
    ('M', 'Y'): -1, ('M', 'V'): 1,
    ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3,
    ('F', 'V'): -1,
    ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
    ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
    ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
    ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
    ('Y', 'Y'): 7, ('Y', 'V'): -1,
    ('V', 'V'): 4
}

# Function to calculate sequence similarity score
def calculate_blosum62_similarity(seq1, seq2, blosum62, gap_penalty=-1, special_char_penalty=-1):
    """
    Calculate the similarity score between two sequences using the BLOSUM62 matrix.
    Sequences are first aligned using edlib.

    Parameters:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        blosum62 (dict): BLOSUM62 matrix as a dictionary.
        gap_penalty (int): Penalty for gaps introduced during alignment.
        special_char_penalty (int): Penalty for special characters like '*'.

    Returns:
        int: Similarity score.
    """
    # Perform global alignment using edlib
    alignment = edlib.align(seq1, seq2, mode="HW", task="path")
    
    # result = edlib.align(rbd_template, sequence, mode="HW", task="path")
    nice = edlib.getNiceAlignment(alignment, seq1, seq2)

    aligned_seq1 = nice['query_aligned']
    aligned_seq2 = nice['target_aligned']
    
    # Calculate similarity score
    score = 0
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-" or b == "-":
            score += gap_penalty  # Penalize gaps
        elif a == "*" or b == "*":
            score += special_char_penalty  # Penalize special characters
        else:
            score += blosum62.get((a, b), blosum62.get((b, a), gap_penalty))  # BLOSUM62 lookup

    return score

# Function to calculate sequence similarity score
def calculate_identity(seq1, seq2):
    """
    Calculate the similarity score between two sequences using the BLOSUM62 matrix.
    Sequences are first aligned using edlib.

    Parameters:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        blosum62 (dict): BLOSUM62 matrix as a dictionary.
        gap_penalty (int): Penalty for gaps introduced during alignment.
        special_char_penalty (int): Penalty for special characters like '*'.

    Returns:
        int: Similarity score.
    """
    # Perform global alignment using edlib
    alignment = edlib.align(seq1, seq2, mode="HW", task="path")
    
    # result = edlib.align(rbd_template, sequence, mode="HW", task="path")
    nice = edlib.getNiceAlignment(alignment, seq1, seq2)

    identity = 0
    aligned_seq1 = nice['query_aligned']
    aligned_seq2 = nice['target_aligned']
    for i in range(len(aligned_seq1)):
        identity += (1 if aligned_seq1[i] == aligned_seq2[i] else 0)
    
    return round(identity/len(aligned_seq1), 2)

def calculate_similarity_matrix(sequences, matrix):
    n = len(sequences)
    similarity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            # if i <= j:  # Only calculate upper triangle
            similarity_matrix[i, j] = 100*calculate_identity(sequences[i], sequences[j])
            # similarity_matrix[i, j] = calculate_blosum62_similarity(sequences[i], sequences[j], BLOSUM62)
            similarity_matrix[j, i] = similarity_matrix[i, j]  # Symmetry

    return similarity_matrix

# Filter missing CDR3
df = df[df['cdr3_aa'].str.len() > 0]

# Cleanup V & J families
df["v_family"] = df["v_call"].str.split("*").str[0]
df["j_family"] = df["j_call"].str.split("*").str[0]

# Create "subject" column
df['subject_id'] = df['sequence_id'].str.split('-', n=1).str[0]
df['subject_id'] = df['subject_id'].str.replace('SBJ', '', regex=False)

# Create an antibody_id
df['antibody_id'] = df['sequence_id'].str.extract('-(.*?)(?=_)', expand=False)

# Create an antibody_name
df['antibody_id'] = df['subject_id'] + df['antibody_id']

# Assign "cluster" column based on the "name" column and dictionary
df['cluster'] = df['antibody_id'].map(cluster_map).fillna(-1)

# Filter out antibodies without a cluster
df = df[ df['cluster'] >= 0]
df.reset_index(drop=True, inplace=True)
print(df.columns)
# Assign chain
df['chain'] = df['locus'].apply(lambda x: 'heavy' if x == 'IGH' else 'light')

# Select a subset of columns
column_list = [
    'subject_id',
    'antibody_id',
    'cluster',
    'chain',
    'sequence',
    'sequence_aa',
    'v_family',
    'j_family',
    'v_identity',
    'd_identity',
    'j_identity'
]
df = df[column_list]
print(df)

# Separate heavy and light chains
df_heavy = df[df['chain'] == 'heavy']
df_light = df[df['chain'] == 'light']

# # Count distributions for V and J genes
# heavy_v_counts = heavy_df['v_family'].value_counts().sort_values()
# light_v_counts = light_df['v_family'].value_counts().sort_values()
# heavy_j_counts = heavy_df['j_family'].value_counts().sort_values()
# light_j_counts = light_df['j_family'].value_counts().sort_values()

# # Plotting
fig, axes = plt.subplots(4, 2, figsize=(8, 11), gridspec_kw={'height_ratios': [0.7, 0.7, 0.7, 1]})
axes = axes.ravel()

# # Top row: Heavy and Light V distributions
# axes[0].barh(heavy_v_counts.index, heavy_v_counts.values, color=color_map['heavy'])
# # axes[0].set_title('Heavy V Germline Distribution')
# axes[0].set_xlabel('Count')
# axes[0].set_ylabel('Heavy V Germline')
# axes[0].text(-0.2, 1.05, 'A', transform=axes[0].transAxes,size=12, weight='bold')

# axes[1].barh(light_v_counts.index, light_v_counts.values, color=color_map['light'])
# # axes[1].set_title('Light V Germline Distribution')
# axes[1].set_xlabel('Count')
# axes[1].set_ylabel('Light V Germline')
# # axes[1].text(-0.2, 1.05, 'B', transform=axes[1].transAxes,size=12, weight='bold')

# # Middle row: Heavy and Light J distributions
# axes[2].barh(heavy_j_counts.index, heavy_j_counts.values, color=color_map['heavy'])
# # axes[2].set_title('Heavy J Germline Distribution')
# axes[2].set_xlabel('Count')
# axes[2].set_ylabel('Heavy J Germline')
# axes[2].text(-0.2, 1.05, 'B', transform=axes[2].transAxes,size=12, weight='bold')

# axes[3].barh(light_j_counts.index, light_j_counts.values, color=color_map['light'])
# # axes[3].set_title('Light J Germline Distribution')
# axes[3].set_xlabel('Count')
# axes[3].set_ylabel('Light J Germline')
# # axes[3].text(-0.2, 1.05, 'D', transform=axes[3].transAxes,size=12, weight='bold')





# Get unique values for "v_family" and "j_family"

# Group by chain, v_gene, and j_gene and get the frequency of combinations
df_heatmap = df_heavy.groupby(['v_family', 'j_family', 'cluster']).size().reset_index(name='frequency')
df_heatmap['percentage'] = df_heatmap['frequency'].transform(lambda x: x / x.sum() * 100)
print(df_heatmap)

v_family_categories = df_heatmap['v_family'].unique()
j_family_categories = df_heatmap['j_family'].unique()

# Map categorical values to integer indices for plotting
df_heatmap['v_family_index'] = df_heatmap['v_family'].map({v: i for i, v in enumerate(v_family_categories)})
df_heatmap['j_family_index'] = df_heatmap['j_family'].map({j: i for i, j in enumerate(j_family_categories)})

# Assign colors for clusters (using a colormap)
cluster_colors = {cluster: color for cluster, color in zip(df_heatmap['cluster'].unique(), plt.cm.get_cmap('viridis', len(df_heatmap['cluster'].unique())).colors)}
print(df_heatmap)
# Create a scatter plot (bubble plot)
axes[0].scatter(df_heatmap['v_family_index'], df_heatmap['j_family_index'], s=df_heatmap['frequency']*10, 
                      c=df_heatmap['cluster'].map(cluster_colors), alpha=0.6)

# Set the axis labels and ticks
axes[0].set_xlabel('V Family')
axes[0].set_ylabel('J Family')
# axis[0].set_title('Bubble Plot of V Family vs J Family')

# Customize x and y axis labels with the category names
axes[0].set_xticks(np.arange(len(v_family_categories)), v_family_categories)
axes[0].set_yticks(np.arange(len(j_family_categories)), j_family_categories)






# Generate the data for the violing plot
# heavy_df.head()
df_violin = df[['chain', 'v_identity','d_identity','j_identity']].melt(id_vars='chain', var_name='germline', value_name='mutations')
df_violin['germline'] = df_violin['germline'].str.replace('_identity', '').str.upper()

# heavy chain violin plot
sns.violinplot(
    data=df_violin[df_violin['chain'] == 'heavy'],
    x='germline',
    y='mutations',
    # palette='muted',
    # inner='box',
    linewidth=0.5,
    inner=None,
    palette=color_genes,
    hue="germline",
    # cut=0,
    ax=axes[4]
)
sns.stripplot(
    data=df_violin[df_violin['chain'] == 'heavy'],
    x='germline',
    y='mutations',
    jitter=True, 
    linewidth=0.5,
    color='#545454',
    alpha=0.6,
    ax=axes[4]
)
axes[4].set_ylabel("Identity (%)")
axes[4].set_xlabel("Gene")
axes[4].set_ylim(60, 105)
axes[4].text(-0.2, 1.05, 'C', transform=axes[4].transAxes,size=12, weight='bold')

# light chain violin plot
sns.violinplot(
    data=df_violin[(df_violin['chain'] == 'light') & (df_violin['germline'] != 'D')],
    x='germline',
    y='mutations',
    palette=color_genes,
    hue="germline",
    inner=None,
    linewidth=0.5,
    ax=axes[5]
)
sns.stripplot(
    data=df_violin[(df_violin['chain'] == 'light') & (df_violin['germline'] != 'D')],
    x='germline',
    y='mutations',
    jitter=True, 
    linewidth=0.5,
    color='#545454',
    alpha=0.6,
    ax=axes[5]
)

axes[5].set_ylabel("Identity (%)")
axes[5].set_xlabel("Gene")
axes[5].set_ylim(60, 105)
# axes[5].text(-0.2, 1.05, 'F', transform=axes[5].transAxes,size=12, weight='bold')

# Group by chain, v_gene, and j_gene and get the frequency of combinations
df_heatmap = df.groupby(['chain', 'v_family', 'j_family']).size().reset_index(name='frequency')
df_heatmap['percentage'] = df_heatmap.groupby('chain')['frequency'].transform(lambda x: x / x.sum() * 100)

axes[6].text(-0.2, 1.05, 'D', transform=axes[6].transAxes,size=12, weight='bold')
for idx, chain in enumerate(['heavy', 'light']):
    # Filter data for the specific chain
    df_chain = df_heatmap[df_heatmap['chain'] == chain]
    
    # Pivot the table for the current chain
    df_chain = df_chain.pivot_table(index='v_family', columns='j_family', values='percentage', fill_value=0)
    
    # Plot the heatmap
    sns.heatmap(df_chain, annot=True, cmap='Blues', ax=axes[6+idx], cbar=False, fmt='.1f')
    axes[6+idx].set_xlabel('J Gene')
    axes[6+idx].set_ylabel('V Gene')

# Adjust layout
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'germline.png'), dpi=600)
plt.close(fig)

# Perform hierarchical clustering
def hierarchical_clustering(similarity_matrix):
    distance_matrix = 1 - similarity_matrix / np.max(similarity_matrix)  # Convert to distances
    np.fill_diagonal(distance_matrix, 0)  # Ensure diagonal is zero for proper clustering

    linkage_matrix = linkage(squareform(distance_matrix), method='average')
    return linkage_matrix

# Create the heatmap figure
fig, axis = plt.subplots(2, 1, figsize=(5, 10))
axis = axis.ravel()

# Calculate sequence similarity matrices
similarity_column = 'sequence_aa'
heavy_similarity = calculate_similarity_matrix(df_heavy[similarity_column].tolist(), BLOSUM62)
light_similarity = calculate_similarity_matrix(df_light[similarity_column].tolist(), BLOSUM62)
heavy_linkage = hierarchical_clustering(heavy_similarity)
light_linkage = hierarchical_clustering(light_similarity)

# Perform hierarchical clustering on the rows and columns
heavy_labels = df_heavy['antibody_id'].tolist()
light_labels = df_light['antibody_id'].tolist()
heavy_row_clusters = linkage(heavy_similarity, method='ward')
heavy_col_clusters = linkage(heavy_similarity.T, method='ward')
light_row_clusters = linkage(light_similarity, method='ward')
light_col_clusters = linkage(light_similarity.T, method='ward')

# Get the ordering of rows and columns based on clustering
heavy_row_order = leaves_list(heavy_row_clusters)
heavy_col_order = leaves_list(heavy_col_clusters)
light_row_order = leaves_list(light_row_clusters)
light_col_order = leaves_list(light_col_clusters)

# Reorder the data according to clustering
clustered_heavy = heavy_similarity[heavy_row_order, :][:, heavy_col_order]
clustered_light = light_similarity[light_row_order, :][:, light_col_order]

# Reorder the labels accordingly
heavy_x_labels = [heavy_labels[i] for i in heavy_col_order]
heavy_y_labels = [heavy_labels[i] for i in heavy_row_order]
light_x_labels = [light_labels[i] for i in light_col_order]
light_y_labels = [light_labels[i] for i in light_row_order]

# Bottom row: Heatmaps with dendrograms
sns.heatmap(
    clustered_heavy, 
    ax=axis[0], 
    cmap='coolwarm', 
    square=True, 
    linewidths=0.1,
    cbar=True, 
    xticklabels=heavy_x_labels, 
    yticklabels=heavy_y_labels,
    cbar_kws=dict(use_gridspec=True,location="right",pad=0.01,shrink=0.83)
)
# axis[0].set_title('Heavy Chain Sequence Similarity')
axis[0].set_xlabel('Antibodies')
axis[0].set_ylabel('Antibodies')

# Reduce label size
axis[0].set_xticklabels(axis[0].get_xticklabels(), fontsize=6)
axis[0].set_yticklabels(axis[0].get_yticklabels(), fontsize=6)
axis[0].text(-0.2, 1.05, 'A', transform=axis[0].transAxes,size=12, weight='bold')

sns.heatmap(
    clustered_light, 
    ax=axis[1], 
    cmap='coolwarm', 
    square=True, 
    linewidths=0.1,
    cbar=True, 
    xticklabels=light_x_labels, 
    yticklabels=light_y_labels,
    cbar_kws=dict(use_gridspec=True,location="right",pad=0.01,shrink=0.83)
)

# axis[0].set_title('Light Chain Sequence Similarity')
axis[1].set_xlabel('Antibodies')
axis[1].set_ylabel('Antibodies')
axis[1].set_xticklabels(axis[1].get_xticklabels(), fontsize=6)
axis[1].set_yticklabels(axis[1].get_yticklabels(), fontsize=6)
axis[1].text(-0.2, 1.05, 'B', transform=axis[1].transAxes,size=12, weight='bold')

# Adjust layout
plt.tight_layout()
plt.savefig(os.path.join(output_path, 'heatmap.png'), dpi=600)
plt.close(fig)