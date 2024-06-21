# Anti-capsule_human_monoclonal_antibodies_protect_against_pandrug-resistant_Klebsiella_pneumoniae
Scripts for the manuscript Anti-capsule human monoclonal antibodies protect against hypervirulent and pandrug-resistant Klebsiella pneumoniae

## Analysis of KL64 and ST147 Kp global diffusion

Download the *"Klebsiella_pneumoniae__kleborate.csv"* and *"Klebsiella_pneumoniae__metadata.csv"* files from [Pathogenwatch](https://pathogen.watch/) following the Public data downloads guidelines at [here](https://cgps.gitbook.io/pathogenwatch/public-data-downloads).

Create and activate the conda environment

```
conda create -p ./env_conda python=3.10 pip ipython
conda activate ./env_conda/
```

Install the required packages

```
conda install pandas
conda install matplotlib
conda install seaborn
conda install geopandas
```

Run the ST_KL_distribution.py script to analyze the data and obtain the figures for 
- the distribution of the STs by year;
- the distribution of the K_locus by year;
- Virulence and resistance score distribution by year.

```
python ST_KL_distribution.py
```

Run the word_map_plot.py to analize the data and obtain the word map figures showing the K locus KL64 geografical distribution for the three time frames till-2010, till-2016 and till-2022.

```
python word_map_plot.py
```
