# Hetairos - Prediction tool for conserved chromosomal clusters of Pfam domain coding regions in microbial genomes
## General info
This repository includes a set of Python 3.6 scripts that were used for the prediction of chromosomal clusters of Pfam domain regions from Bacterial and Archaeal genomes.

The name hetairos comes from the ancient greek word ἑταῖρος (hetaîros, genitive ἑταῖρου), a noun meaning comrade, companion, partner, friend; pupil, disciple; member of a religious guild; (rarely) lover, member of a pair of lovers; (in the plural) the guards i.e. the cavalry of the Macedonian kings; (as an adjective) associate of. Just like two Pfam domains that cluster throughout evolution. :) 
	
## Setup
The scripts were written to be used with Python v3.6. The required libraries are listed in the file requirements.txt.
Installation instructions are provided for anaconda. From the command line:

```
$ git clone https://github.com/selenocysteine/conserved-chromosomal-clusters.git
$ cd conserved-chromosomal-clusters
$ conda create -n hetairos python=3.6
$ conda activate hetairos
$ pip install -r requirements.txt
```

## Running the code
The input files required for running this analysis are .gbk files from each microbial species, coupled with Pfam annotations files as returned from HMMER 3.0. 
The analysis can be run from src/main.py with the following parameters:

```
usage: main.py [-h] [-t THREADS] [-i DATA_DIR] [-o OUT_DIR] [--fit]
               [--predict] [--weight] [--subclades] [--permutations]
               [--nperm NPERM] [--plots] [--evalue EVALUE] [--phi PHI_0]
               [--lambda LAMBD_0] [--bootstrap] [--sensitivity]
               [--tree TREE_FILE] [--save_pairs] [--non_go]
               [--pfam PFAM_ANNOTATIONS_FILENAME_SUFFIX] [--force] [--verbose]

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS            Number of threads to use (default=1).
  -i DATA_DIR           Input files directory.
  -o OUT_DIR            Output files directory.
  --fit                 Estimate new clustering distances parameters using the
                        organisms in the input mapping_tables directory
                        (optional)
  --predict             Predict clustering probabilities for all the pairs of
                        Pfam domains from the organisms in the input
                        mapping_tables directory (default=False).
  --weight              While computing conserved clustering probabilities,
                        weight each species accordingto its position on a
                        phylogenetic tree (default=False)
  --subclades           Compute a set of subclade-specific clustering
                        probabilities instead of just a single probability for
                        the whole dataset (default=False)
  --permutations        Compute random permutations of gene distances to
                        assess significances
  --nperm NPERM         Number of permutations (default=1).
  --plots               Generate plots for the fitting step (default=False).
  --evalue EVALUE       E-value threshold for Pfam predict (default=0.001).
  --phi PHI_0           Initial value for phi (default=0.001).
  --lambda LAMBD_0      Initial value for the rate parameter (default=uses a
                        moment based estimator)
  --bootstrap           Perform bootstrap stability analysis of the parameters
                        (default=False)
  --sensitivity         Perform sensitivity analysis with mean parameters
                        (default=False)
  --tree TREE_FILE      Path to the file containing the phylogenetic tree of
                        the analysed species
  --save_pairs          Save organism-specific domain pair distances values on
                        disk (default=False)
  --non_go              During the fit, compute pair distances of non-GO
                        related domains (default=False)
  --pfam PFAM_ANNOTATIONS_FILENAME_SUFFIX
                        Filename suffix for the Pfam input files
                        (default=pfam_annotations)
  --force               Completely rewrite eventual previously existing input
                        files (default=False)
  --verbose             Print detailed comments (default=False)

```
