import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import geopandas
import numpy as np

## --- Parameters --- ##
replace_unknown_kl64 = False
point_color = 'red'
point_size = 18
alpha_jitter = 1. # .1
markersize=point_size
edgecolor='black'
linewidth=0.2
world_background_color = sns.color_palette("pastel")[5]
# world_background_color = sns.color_palette("hls", 8)[4]
# world_background_color = sns.color_palette("pastel")[7]
# world_background_color = sns.color_palette("pastel")[2]
## --- [END] Parameters --- ##

## Data loading and merging
file_path = "/data"

# load kleborate file
kleborate_path = os.path.join(file_path, 'Klebsiella_pneumoniae__kleborate.csv')
kleborate_df = pd.read_csv(kleborate_path, sep=",")

kleborate_var_list = ['Genome ID',
 'Genome Name',
 'strain',
 'species',
 'ST',
 'virulence_score',
 'resistance_score',
 'num_resistance_classes',
 'num_resistance_genes',
 'K_locus',
 'K_type',
 'O_locus',
 'O_type']

# select the columns of interest
kleborate_df = kleborate_df[kleborate_var_list]

# # -- K locus KL64 -- #
print("number of KL64: ", kleborate_df['K_locus'].value_counts()) # number of KL64: 3742.

# load metadata file
metadata_path = os.path.join(file_path, 'Klebsiella_pneumoniae__metadata.csv')
metadata_df = pd.read_csv(metadata_path, sep=",")

metadata_var_list = ['id',
 'displayname',
 'latitude',
 'longitude',
 'year',
 'Host',
 'Collection date',
 'City/region',
 'Country']

# select the columns of interest
metadata_df = metadata_df[metadata_var_list]

# merge kleborate and metadata files
df_merged = pd.merge(kleborate_df, metadata_df, left_on='Genome ID', right_on='id')

# drop the columns that are not needed
df_merged = df_merged.drop(columns=['id', 'displayname'])

# drop the rows with missing values in year and ST columns
df_merged = df_merged.dropna(subset=['year', 'ST', 'K_locus'])

# cast the year column to integer
df_merged['year'] = df_merged['year'].astype(int)

# # -- Remove this filter to compatibility with Plot sent to the authors with AK email the 21May24 -- #
# # filter the rows with year greater than 2010
# df_merged = df_merged[df_merged['year'] >= 2010]

# remove the rows with year 2023 and 2024
df_merged = df_merged[df_merged['year'] <= 2022]
print(df_merged.shape) # (39536, 20)


# add a bit of random jitter to the latitude and longitude columns to avoid overlapping of the points
sigma = alpha_jitter
df_merged['latitude'] = df_merged['latitude'].apply(lambda x: x + np.random.normal(0, sigma))
df_merged['longitude'] = df_merged['longitude'].apply(lambda x: x + np.random.normal(0, sigma))

# load the world map
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))


# -- Create a word map plot with K locus KL64 for the year 2010 -- #
# filter the rows with year <= 2010. Cumulative number of isolates till 2010
df_merged_2010 = df_merged[df_merged['year'] <= 2010]
# filter the rows with K_locus KL64
df_merged_2010_kl64 = df_merged_2010[df_merged_2010['K_locus'] == 'KL64']

# create a geodataframe from the merged dataframe
gdf_2010_kl64 = geopandas.GeoDataFrame(df_merged_2010_kl64, geometry=geopandas.points_from_xy(df_merged_2010_kl64.longitude, df_merged_2010_kl64.latitude))
# plot the world map
fig, ax = plt.subplots(figsize=(15, 8), dpi=300)
# plot a basic map of the world
world.plot(
    ax=ax,
    color=world_background_color,
    edgecolor="black",
    alpha=0.3
)
gdf_2010_kl64.plot(ax=ax, color=point_color, markersize=point_size, edgecolor=edgecolor, linewidth=linewidth)
# remove the xticks and yticks
plt.xticks([])
plt.yticks([])
# reduce the space between the plot and the edge of the figure
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.title(f'Cumulative number of genomes reported in Pathogen watch for Klebsiella pneumoniae with K64 capsule up to December 2010. \nTotal number of isolates: {len(gdf_2010_kl64)}.', fontsize=12)
plt.savefig(f'world_map_plot_2010_kl64_replace_unknown_kl64{replace_unknown_kl64}.png')
plt.close()
# -- [END] Create a word map plot with K locus KL64 for the year 2010 -- #


# -- Create a word map plot with K locus KL64 for the year 2016 -- #
# filter the rows with year <= 2016
df_merged_2016 = df_merged[df_merged['year'] <= 2016]
# filter the rows with K_locus KL64
df_merged_2016_kl64 = df_merged_2016[df_merged_2016['K_locus'] == 'KL64']

# create a geodataframe from the merged dataframe
gdf_2016_kl64 = geopandas.GeoDataFrame(df_merged_2016_kl64, geometry=geopandas.points_from_xy(df_merged_2016_kl64.longitude, df_merged_2016_kl64.latitude))
# plot the world map
fig, ax = plt.subplots(figsize=(15, 8), dpi=300)
# plot a basic map of the world
world.plot(
    ax=ax,
    color=world_background_color,
    edgecolor="black",
    alpha=0.3
)
gdf_2016_kl64.plot(ax=ax, color=point_color, markersize=point_size, edgecolor=edgecolor, linewidth=linewidth)
# remove the xticks and yticks
plt.xticks([])
plt.yticks([])
# reduce the space between the plot and the edge of the figure
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.title(f'Cumulative number of genomes reported in Pathogen watch for Klebsiella pneumoniae with K64 capsule up to December 2016. \nTotal number of isolates: {len(gdf_2016_kl64)}.', fontsize=12)
plt.savefig(f'world_map_plot_2016_kl64_replace_unknown_kl64{replace_unknown_kl64}.png')
plt.close()
# -- [END] Create a word map plot with K locus KL64 for the year 2016 -- #


# -- Create a word map plot with K locus KL64 for the year 2022 -- #
# filter the rows with year <= 2022
df_merged_2022 = df_merged[df_merged['year'] <= 2022]
# filter the rows with K_locus KL64
df_merged_2022_kl64 = df_merged_2022[df_merged_2022['K_locus'] == 'KL64']

# create a geodataframe from the merged dataframe
gdf_2022_kl64 = geopandas.GeoDataFrame(df_merged_2022_kl64, geometry=geopandas.points_from_xy(df_merged_2022_kl64.longitude, df_merged_2022_kl64.latitude))
# plot the world map
fig, ax = plt.subplots(figsize=(15, 8), dpi=300)
# plot a basic map of the world
world.plot(
    ax=ax,
    color=world_background_color,
    edgecolor="black",
    alpha=0.3
)
gdf_2022_kl64.plot(ax=ax, color=point_color, markersize=point_size, edgecolor=edgecolor, linewidth=linewidth)
# remove the xticks and yticks
plt.xticks([])
plt.yticks([])
# plt.tight_layout()
# reduce the space between the plot and the edge of the figure
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.title(f'Cumulative number of genomes reported in Pathogen watch for Klebsiella pneumoniae with K64 capsule up to December 2022. \nTotal number of isolates: {len(gdf_2022_kl64)}.', fontsize=12)
plt.savefig(f'world_map_plot_2022_kl64_replace_unknown_kl64{replace_unknown_kl64}.png')
plt.close()
# -- [END] Create a word map plot with K locus KL64 for the year 2022 -- #

print('Done!')

