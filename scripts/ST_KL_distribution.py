import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


# ---------------------------------------
# normalize the STs distribution in percentage by year
normalise_precent = True
compute_over_the_most_frequent = False

# Siena strain characterisation
siena_ST = 'ST147'
siena_K_locus = 'KL64'

file_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

# Create the output directory
result_path = os.path.join(file_path, 'results')
if not os.path.exists(result_path):
    os.mkdir(result_path)
    
# ---------------------------------------

## Data loading and merging
# load kleborate file
kleborate_path = os.path.join(file_path, 'data', 'klebsiella_pneumoniae__kleborate.csv')
kleborate_df = pd.read_csv(kleborate_path, sep=",")

kleborate_var_list = [
    'Genome ID',
    'Genome Name',
    #  'strain',
    'species',
    'ST',
    'virulence_score',
    'resistance_score',
    'resistance_class_count',
    'resistance_gene_count',
    'K_locus',
    'K_type',
    'O_locus',
    'O_type'
]

# select the columns of interest
kleborate_df = kleborate_df[kleborate_var_list]

# # -- K locus KL64 -- #
print("number of KL64: ", kleborate_df['K_locus'].value_counts()) # number of KL64: 3742.


# load metadata file
metadata_path = os.path.join(file_path, 'data', 'klebsiella_pneumoniae__metadata.csv')
metadata_df = pd.read_csv(metadata_path, sep=",")

metadata_var_list = ['id',
 'displayname',
 'latitude',
 'longitude',
 'year',
 'Host',
 'Collection date',
 'City/region',
 'Host',
 'Environmental sample',
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

# filter the rows with year greater than 2010
df_merged = df_merged[df_merged['year'] >= 2010]

# remove the rows with year 2023 and 2024
df_merged = df_merged[df_merged['year'] <= 2022]

## Prevalence information
# List of the most frequent STs
print(df_merged['ST'].value_counts().head(10))
most_frequent_STs = df_merged['ST'].value_counts().head(10).index.tolist()

# List of the most frequent K_locus
print(df_merged['K_locus'].value_counts().head(11))
most_frequent_K_locus = df_merged['K_locus'].value_counts().head(11).index.tolist()
most_frequent_K_locus = [_ for _ in most_frequent_K_locus if 'unknown' not in _]

# Make sure some of the most relevant are present & remove duplicates
most_frequent_K_locus.extend(['KL1', 'KL2', 'KL64'])
most_frequent_K_locus = list(set(most_frequent_K_locus))

# List of the most frequent O_locus
print(df_merged['O_locus'].value_counts().head(10))


# get seaborn default color palette
palette = sns.color_palette()

# --------------------------------------- Plot of the distribution of the STs by year

st_by_year_df = df_merged.groupby(by=['year', 'ST']).size().unstack().reset_index()#.drop(columns='0')
st_by_year_df = st_by_year_df.fillna(0)


if compute_over_the_most_frequent:
    # Select the most frequent STs
    st_by_year_df = st_by_year_df[['year'] + most_frequent_STs]
    if normalise_precent:
        st_by_year_df = st_by_year_df.set_index('year')
        st_by_year_df = st_by_year_df.div(st_by_year_df.sum(axis=1), axis=0) * 100
        st_by_year_df = st_by_year_df.reset_index()
else:    
    if normalise_precent:
        st_by_year_df = st_by_year_df.set_index('year')
        st_by_year_df = st_by_year_df.div(st_by_year_df.sum(axis=1), axis=0) * 100
        st_by_year_df = st_by_year_df.reset_index()
    # Select the most frequent STs
    st_by_year_df = st_by_year_df[['year'] + most_frequent_STs]

# get the siena_ST distribution by year
siena_st_df = st_by_year_df[[_ for _ in list(st_by_year_df) if _  in ['year',siena_ST]]]
most_frequent_st_df = st_by_year_df[[_ for _ in list(st_by_year_df) if _  != siena_ST]]

# --------------------------------------- Draw: Yearly distribution of the percentage of isolates for the ten most frequent STs
# Figure
plt.figure(dpi=300)
sns.set(style="whitegrid") #{'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
ax = sns.lineplot(x='year', y='value', hue='ST', data=pd.melt(siena_st_df, ['year']), legend='brief',lw=1.5)
ax1 = sns.lineplot(x='year', y='value', hue='ST', data=pd.melt(most_frequent_st_df, ['year']), legend='brief',lw=.8, alpha=.7, palette=palette[1:], linestyle='dashdot')
plt.xlabel('Year')
# reduce the legend font to 'x-small'
sns.move_legend(ax, "upper left", fontsize=8, title=None, ncol=2)
sns.move_legend(ax1, "upper left", fontsize=8, title=None, ncol=2)

if normalise_precent:
    plt.ylabel('Isolates (%)')
    plt.title('Percentage of isolates by year')
    if compute_over_the_most_frequent:
        plt.title('Yearly distribution of the percentage of isolates for the ten most frequent STs\n(percentage computed over the ten most frequent STs)', fontsize=11)
        plt.savefig(os.path.join(result_path, 'number_of_isolates_by_year_and_ST_lineplot_percentage_over_all_the_STs_over_most_frequent.png'))
    else:
        plt.title('Yearly distribution of the percentage of isolates for the ten most frequent STs\n(percentage computed over all the STs)', fontsize=11)
        plt.title(" ")
        plt.savefig(os.path.join(result_path, 'number_of_isolates_by_year_and_ST_lineplot_percentage_over_all_the_STs.png'))
else:
    plt.ylabel('Number of isolates')
    plt.title('Number of isolates by year')
    plt.savefig(os.path.join(result_path, 'number_of_isolates_by_year_and_ST_lineplot.png'))
plt.close()
# ---------------------------------------

# --------------------------------------- Plot of the distribution of the K_locus by year
# Get the K_locus distribution by year
k_locus_by_year_df = df_merged.groupby(by=['year', 'K_locus']).size().unstack().reset_index()
k_locus_by_year_df = k_locus_by_year_df.fillna(0)

if compute_over_the_most_frequent:
    # Select the most frequent K_locus
    k_locus_by_year_df = k_locus_by_year_df[['year'] + most_frequent_K_locus]
    if normalise_precent:
        k_locus_by_year_df = k_locus_by_year_df.set_index('year')
        k_locus_by_year_df = k_locus_by_year_df.div(k_locus_by_year_df.sum(axis=1), axis=0) * 100
        k_locus_by_year_df = k_locus_by_year_df.reset_index()
else:
    # normalize the K_locus distribution in percentage by year
    if normalise_precent:
        k_locus_by_year_df = k_locus_by_year_df.set_index('year')
        k_locus_by_year_df = k_locus_by_year_df.div(k_locus_by_year_df.sum(axis=1), axis=0) * 100
        k_locus_by_year_df = k_locus_by_year_df.reset_index()
    # Select the most frequent K_locus
    k_locus_by_year_df = k_locus_by_year_df[['year'] + most_frequent_K_locus]

# get the siena_K_locus distribution by year
siena_k_locus_df = k_locus_by_year_df[[_ for _ in list(k_locus_by_year_df) if _  in ['year',siena_K_locus]]]
most_frequent_k_locus_df = k_locus_by_year_df[[_ for _ in list(k_locus_by_year_df) if _  != siena_K_locus]]


# --------------------------------------- Draw: Yearly distribution of the percentage of isolates for the ten most frequent K_locus
# Figure
plt.figure(dpi=300)
sns.set(style="whitegrid") #{'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
ax = sns.lineplot(x='year', y='value', hue='K_locus', data=pd.melt(siena_k_locus_df, ['year']), legend='brief',lw=1.5)
ax1 = sns.lineplot(x='year', y='value', hue='K_locus', data=pd.melt(most_frequent_k_locus_df, ['year']), legend='brief',lw=.8, alpha=.7, palette=palette[1:], linestyle='dashdot')
plt.xlabel('Year')

# reduce the legend font to 'x-small'
sns.move_legend(ax, "upper left", fontsize=8, title=None, ncol=2)
sns.move_legend(ax1, "upper left", fontsize=8, title=None, ncol=2)
if normalise_precent:
    plt.ylabel('Isolates (%)')
    if compute_over_the_most_frequent:
        plt.title('Yearly distribution of the percentage of isolates for the ten most frequent K_locus\n(percentage computed over the ten most frequent K_locus)', fontsize=11)
        plt.savefig(os.path.join(result_path, 'number_of_isolates_by_year_and_K_locus_lineplot_percentage_over_most_frequent.png'))

    else:
        plt.title('Yearly distribution of the percentage of isolates for the ten most frequent K_locus\n(percentage computed over all the K_locus)', fontsize=11)
        plt.title(" ")
        plt.savefig(os.path.join(result_path, 'number_of_isolates_by_year_and_K_locus_lineplot_percentage_over_all_the_STs.png'))
else:
    plt.ylabel('Number of isolates')
    plt.title('Number of isolates by year')
    plt.savefig(os.path.join(result_path, 'number_of_isolates_by_year_and_K_locus_lineplot.png'))

plt.close()
# ---------------------------------------

# --------------------------------------- Virulence and resistance score distribution by year
# group the strains by year and ST and counting the number of item with virulence score > 3 and resistance scores > 3
df_merged['virulence_score_gt_3'] = df_merged['virulence_score'] >= 3
df_merged['resistance_score_gt_2'] = df_merged['resistance_score'] >= 2

# group the strains by year and ST and counting the number of item with virulence score < 3 and resistance scores < 3
df_merged['virulence_score_ls_3'] = df_merged['virulence_score'] < 3
df_merged['resistance_score_ls_2'] = df_merged['resistance_score'] < 2

# normalize the virulence and resistance scores ge 3 in percentage with rispect to the number of resistance scores ls 3 by year
if normalise_precent:
    virulence_resistance_df = df_merged.groupby(['year', 'ST'])[['virulence_score_gt_3', 'virulence_score_ls_3', 'resistance_score_gt_2', 'resistance_score_ls_2']].sum().reset_index()
    virulence_resistance_df = virulence_resistance_df.set_index('year')
    # normalize the virulence and resistance scores ge 3 in percentage with rispect to the number of resistance scores ls 3 by year
    virulence_resistance_df['virulence_score_gt_3'] = virulence_resistance_df['virulence_score_gt_3'] / (virulence_resistance_df['virulence_score_gt_3'] + virulence_resistance_df['virulence_score_ls_3']) * 100
    virulence_resistance_df['resistance_score_gt_2'] = virulence_resistance_df['resistance_score_gt_2'] / (virulence_resistance_df['resistance_score_gt_2'] + virulence_resistance_df['resistance_score_ls_2']) * 100
    virulence_resistance_df = virulence_resistance_df.reset_index()
else:
    # get the number of isolates with virulence score > 2 and resistance score > 3 by year and ST
    virulence_resistance_df = df_merged.groupby(['year', 'ST'])[['virulence_score_gt_3', 'resistance_score_gt_2']].sum().reset_index()

# select the siena_ST and the most frequent STs
virulence_resistance_df_siena = virulence_resistance_df[virulence_resistance_df['ST'] == siena_ST]
virulence_resistance_df_most_frequent = virulence_resistance_df[virulence_resistance_df['ST'].map(lambda x:x in most_frequent_STs)]
virulence_resistance_df_most_frequent = virulence_resistance_df_most_frequent[virulence_resistance_df_most_frequent['ST'] != siena_ST]

# sort the virulence_resistance_df_most_frequent by the exact same order of the element in most_frequent_STs list
virulence_resistance_df_most_frequent['custom_order'] = virulence_resistance_df_most_frequent['ST'].map(lambda x: most_frequent_STs.index(x))
virulence_resistance_df_most_frequent = virulence_resistance_df_most_frequent.sort_values(by='custom_order')


# --------------------------------------- Draw: Yearly distribution of the percentage of isolates with virulence score > 3 and resistance score > 3
# --------------------------------------- Draw the Figure
# Figure
plt.figure(dpi=300)
sns.set(style="whitegrid") #{'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
ax = sns.lineplot(x='year', y='virulence_score_gt_3', hue='ST', data=virulence_resistance_df_siena, legend='brief', lw=1.5)
ax1 = sns.lineplot(x='year', y='virulence_score_gt_3', hue='ST', data=virulence_resistance_df_most_frequent, legend='brief', lw=.8, alpha=.7, palette=palette[1:], linestyle='dashdot')
plt.xlabel('Year')
sns.move_legend(ax, "upper left", fontsize=8, title=None, ncol=2)
sns.move_legend(ax1, "upper left", fontsize=8, title=None, ncol=2)
if normalise_precent:
    plt.ylabel('Isolates (%)')
    plt.title('Yearly distribution of the percentage of isolates with\n Virulence Score greater than 2 for the ten most frequent STs\n(percentage computed over all the STs)', fontsize=11)
    plt.title(" ")
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_virulence_score_gt_3_by_year_percentage_over_all_the_STs.png'))
else:
    plt.ylabel('Number of isolates')
    plt.title('Yearly distribution of the number of isolates with\n Virulence Score greater than 2 for the ten most frequent STs')
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_virulence_score_gt_3_by_year.png'))
plt.close()
# ---------------------------------------

# --------------------------------------- Draw the Figure
# Figure
plt.figure(dpi=300)
sns.set(style="whitegrid") #{'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
ax = sns.lineplot(x='year', y='resistance_score_gt_2', hue='ST', data=virulence_resistance_df_siena, legend='brief', lw=1.5)
ax1 = sns.lineplot(x='year', y='resistance_score_gt_2', hue='ST', data=virulence_resistance_df_most_frequent, legend='brief', lw=.8, alpha=.7, palette=palette[1:], linestyle='dashdot')
plt.xlabel('Year')
sns.move_legend(ax, "upper left", fontsize=8, title=None, ncol=2)
sns.move_legend(ax1, "upper left", fontsize=8, title=None, ncol=2)
if normalise_precent:
    plt.ylabel('Isolates (%)')
    plt.title('Yearly distribution of the percentage of isolates with\n Resistance Score 2-3 for the ten most frequent ST', fontsize=11)
    plt.title(" ")
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_resistance_score_gt_2_by_year_percentage_over_all_the_STs.png'))
else:
    plt.ylabel('Number of isolates')
    plt.title('Yearly distribution of the number of isolates with\n Resistance Score greater than 2 for the ten most frequent STs')
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_resistance_score_gt_2_by_year.png'))
plt.close()
# ---------------------------------------

# --------------------------------------- K_locus distribution by year
if normalise_precent:
    virulence_resistance_df = df_merged.groupby(['year', 'K_locus'])[['virulence_score_gt_3', 'virulence_score_ls_3', 'resistance_score_gt_2', 'resistance_score_ls_2']].sum().reset_index()
    virulence_resistance_df = virulence_resistance_df.set_index('year')
    # normalize the virulence and resistance scores ge 3 in percentage with rispect to the number of resistance scores ls 3 by year
    virulence_resistance_df['virulence_score_gt_3'] = virulence_resistance_df['virulence_score_gt_3'] / (virulence_resistance_df['virulence_score_gt_3'] + virulence_resistance_df['virulence_score_ls_3']) * 100
    virulence_resistance_df['resistance_score_gt_2'] = virulence_resistance_df['resistance_score_gt_2'] / (virulence_resistance_df['resistance_score_gt_2'] + virulence_resistance_df['resistance_score_ls_2']) * 100
    virulence_resistance_df = virulence_resistance_df.reset_index()
else:
    virulence_resistance_df = df_merged.groupby(['year', 'K_locus'])[['virulence_score_gt_3', 'resistance_score_gt_2']].sum().reset_index()

virulence_resistance_df_siena = virulence_resistance_df[virulence_resistance_df['K_locus'] == siena_K_locus]
virulence_resistance_df_most_frequent = virulence_resistance_df[virulence_resistance_df['K_locus'].map(lambda x:x in most_frequent_K_locus)]
virulence_resistance_df_most_frequent = virulence_resistance_df_most_frequent[virulence_resistance_df_most_frequent['K_locus'] != siena_K_locus]

# sort the virulence_resistance_df_most_frequent by the exact same order of the element in most_frequent_K_locus list
virulence_resistance_df_most_frequent['custom_order'] = virulence_resistance_df_most_frequent['K_locus'].map(lambda x: most_frequent_K_locus.index(x))
virulence_resistance_df_most_frequent = virulence_resistance_df_most_frequent.sort_values(by='custom_order')


# --------------------------------------- Draw the Figure
# Figure
plt.figure(dpi=300)
sns.set(style="whitegrid") #{'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
ax = sns.lineplot(x='year', y='virulence_score_gt_3', hue='K_locus', data=virulence_resistance_df_siena, legend='brief', lw=1.5)
ax1 = sns.lineplot(x='year', y='virulence_score_gt_3', hue='K_locus', data=virulence_resistance_df_most_frequent, legend='brief', lw=.8, alpha=.7, palette=palette[1:], linestyle='dashdot')
plt.xlabel('Year')
sns.move_legend(ax, "upper left", fontsize=8, title=None, ncol=2)
sns.move_legend(ax1, "upper left", fontsize=8, title=None, ncol=2)
if normalise_precent:
    plt.ylabel('Isolates (%)')
    plt.title('Yearly distribution of the percentage of isolates with\n Virulence Score greater than 2 for the ten most frequent K_locus\n(percentage computed over all the K_locus)', fontsize=11)
    plt.title(" ")
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_virulence_score_gt_3_by_year_percentage_KL_over_all_the_STs.png'))
else:
    plt.ylabel('Number of isolates')
    plt.title('Yearly distribution of the number of isolates with\n Virulence Score greater than 2 for the ten most frequent K_locus')
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_virulence_score_gt_3_by_year_KL.png'))
plt.close()
# ---------------------------------------

# --------------------------------------- Draw the Figure
# Figure
plt.figure(dpi=300)
sns.set(style="whitegrid") #{'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
ax = sns.lineplot(x='year', y='resistance_score_gt_2', hue='K_locus', data=virulence_resistance_df_siena, legend='brief', lw=1.5)
ax1 = sns.lineplot(x='year', y='resistance_score_gt_2', hue='K_locus', data=virulence_resistance_df_most_frequent, legend='brief', lw=.8, alpha=.7, palette=palette[1:], linestyle='dashdot')
plt.xlabel('Year')
sns.move_legend(ax, "upper left", fontsize=8, title=None, ncol=2)
sns.move_legend(ax1, "upper left", fontsize=8, title=None, ncol=2)
if normalise_precent:
    plt.ylabel('Isolates (%)')
    plt.title('Yearly distribution of the percentage of isolates with\n Resistance Score 2-3 for the ten most frequent K_locus', fontsize=11)
    plt.title(" ")
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_resistance_score_gt_2_by_year_percentage_KL_over_all_the_STs.png'))
else:
    plt.ylabel('Number of isolates')
    plt.title('Yearly distribution of the number of isolates with\n Resistance Score 2-3 for the ten most frequent K_locus')
    plt.savefig(os.path.join(result_path, 'number_of_isolates_with_resistance_score_gt_2_by_year_KL.png'))
plt.close()

# ---------------------------------------

# Calculate a stacked graph shoing the number of KL64, KL1 and KL2 by year and ST

# df_KL_ST = df_merged.groupby(['year', 'K_locus', 'ST']).size().reset_index(name='count')

# Keep only data for KL64, KL1 and KL2
# df_KL_ST = df_KL_ST[df_KL_ST['K_locus'].isin(['KL64','KL1','KL2'])]

# Select only the specified columns
# df_KL_ST = df_KL_ST[ [ 'year', 'K_locus', 'ST', 'count' ]]

# print(df_KL_ST)


# Count the total number of rows grouped by "year"
total_per_year = df_merged.groupby('year').size().reset_index(name='total_count')

# Step 2: Group by "year", "locus", "st" and count occurrences, then calculate percentages
df_KL_ST = df_merged.groupby(['year', 'K_locus', 'ST']).size().reset_index(name='raw_count')

# Select only the specified columns
df_KL_ST = df_KL_ST[ [ 'year', 'K_locus', 'ST', 'raw_count' ]]

# Merge with total counts to calculate the percentage
df_KL_ST = pd.merge(df_KL_ST, total_per_year, on='year')

# Calculate the percentage
df_KL_ST['percentage'] = (df_KL_ST['raw_count'] / df_KL_ST['total_count']) * 100


# Create a figure with 3 rows and 1 column of subplots
locus_list = ['KL64','KL1','KL2']
ncols = 2
nrows = len(locus_list)
fig, axes = plt.subplots(nrows=len(locus_list), ncols=2, figsize=(ncols*6, nrows*4))
if len(locus_list) == 1:
    axes = [axes]
axes = axes.ravel()

for idx, locus in enumerate(locus_list):
    locus_data = df_KL_ST[df_KL_ST['K_locus'] == locus]
    
    # Group by "ST" and sum the "raw_count"
    grouped_st = locus_data.groupby('ST')['raw_count'].sum()

    # Get the top 10 "ST" by raw_count
    top_10_st = grouped_st.nlargest(10).index.tolist()

    # # Get top N most represented st per year and locus
    # top_n = 10
    # top_st_per_group = locus_data.groupby(['ST'])['ST', 'raw_count'].nlargest(top_n).reset_index(drop=True)
    # print('locus')
    # print(top_st_per_group)
    # sys.exit(1)
    
    # Mark 'st' as 'others' if it's not in the top N
    locus_data['ST'] = locus_data.apply(lambda row: row['ST'] if row['ST'] in top_10_st else 'others', axis=1)

    # Re-group by year, locus, and 'st' (including 'others')
    locus_data = locus_data.groupby(['year', 'K_locus', 'ST'], as_index=False).agg({'raw_count': 'sum'})

    # Recalculate percentages
    locus_data = locus_data.merge(total_per_year, on='year')
    locus_data['percentage'] = (locus_data['raw_count'] / locus_data['total_count']) * 100
    
    # Calculate percentage by year
    locus_by_year = locus_data.groupby('year')['raw_count'].sum().reset_index()
    locus_by_year.rename(columns={'raw_count': 'total_year'}, inplace=True)
    locus_data = pd.merge(locus_data, locus_by_year, on='year', how='left')
    locus_data['percentage_year'] = (locus_data['raw_count'] / locus_data['total_year']) * 100
    
    # locus_data = locus_data.groupby(['year'], as_index=False).agg({'raw_count': 'sum'})
    print(locus_data)
    

    # Pivot the data to have 'st' as columns for stacked bar plot and normalised by total_count
    pivot_data = locus_data.pivot(index='year', columns='ST', values='percentage').fillna(0)

    # Plotting on the specified axis
    axis_idx = idx*2
    pivot_data.plot(kind='bar', stacked=True, ax=axes[axis_idx+0])
    axes[axis_idx+0].set_title(f'{locus} locus')
    axes[axis_idx+0].set_xlabel('Year')
    axes[axis_idx+0].set_ylabel('Isolates (%)')
    axes[axis_idx+0].legend(title=None, ncol=2)
    
    # Plot normalised by year
    pivot_data = locus_data.pivot(index='year', columns='ST', values='percentage_year').fillna(0)
    pivot_data.plot(kind='bar', stacked=True, ax=axes[axis_idx+1], legend=False)
    axes[axis_idx+1].set_title(f'{locus} locus normalised by year')
    axes[axis_idx+1].set_xlabel('Year')
    axes[axis_idx+1].set_ylabel('Isolates (%)')
    # axes[axis_idx+1].legend(title=None, ncol=2)

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the plot to a file
plt.savefig(os.path.join(result_path, 'KL_by_ST_by_year_percentage_KL_over_all_the_STs.png'), dpi=300)


# ---------------------------------------

print('\n\nDone!')
print()
