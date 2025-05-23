# K_pneumoniae-anti_capsule_mabs
Scripts for the manuscript "Anti-capsule human monoclonal antibodies protect against hypervirulent and pandrug-resistant Klebsiella pneumoniae"

## Content ##
```
. 
├── README.md   # This file
├── scripts     # Scripts to generate images
└── data
    ├── blast_ko-locus.tsv      # blastn results for the genes composing K-locus and O-locus
    ├── kleborate.txt           # kleborate annotation results
    ├── core_snp_alignment.tree # Newick tree generated using core genes
    └── antibodies_igblast.tsv   # IgBlast Annotation of the antibodies
```

## Analysis of KL64 and ST147 Kp global diffusion

Run the *download.sh* script that will download data from [Pathogenwatch](https://pathogen.watch/) following the Public data downloads guidelines at [here](https://cgps.gitbook.io/pathogenwatch/public-data-downloads). Files will be placed in the /data folder

Create and activate the python environment

```
python3 -m venv env
source ./env/bin/activate
pip install -r requirements.txt
```

Run the ST_KL_distribution.py script to analyze the data and obtain the figures for 
- the distribution of the STs by year;
- the distribution of the K_locus by year;
- Virulence and resistance score distribution by year.

```
python scripts/ST_KL_distribution.py
```

Run the word_map_plot.py to analize the data and obtain the word map figures showing the K locus KL64 geografical distribution for the three time frames till-2010, till-2016 and till-2022.

```
python scripts/word_map_plot.py
```

## Analysis of Antibody germline

In this script the distribution of heavy and light germline will be calculated and plotted, together with a distance matrix representing the sequence similarity of both heavy and light chains.

    Rscript ./scripts/germline.R -i ./data/antibodies_igblast.tsv -o ./results/germline.pdf

## K-Locus and O-locus analysis

The analysis of the genes included in the K-locus and O-locus of studied strains is performed using `map.R` script. The script assembles the information obtained using Kleborate and blastn:

    Rscript ./scripts/map.R -i ./data/blast_ko-locus.tsv -o ./results/map.png

## Phylogenetic analysis of selected K.pneumoniae strains

Phylogenetic analysis was performed using the pan-genome analysis of the selected strains, and plotted with the following script:

    Rscript ./scripts/tree.R -t ./data/fasttree_cgMLST.nwk -k ./data/kleborate.txt -o ./results/tree.png


