# RPxSuite
RPxSuite is a dedicated, stand-alone, tool for the identification and clustering of marker genes in shotgun metagenomic assemblies.

The publication supporting this tool is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/647511v1)

## Installation with pip
```
$ git clone git@github.com:alexcritschristoph/RPxSuite.git

$ cd RPxSuite

$ pip install .
```

## Quick start
```
$ RPxSuite.py *.fasta -g Ribosomal_L6
```

## Output files

Several files will be produced in the output directory which are described below

#### Finalized files

* `longest.contigs.fasta` has the longest scaffold housing a gene from each gene cluster

* `centroids.fasta` lists the centroid genes in .fasta format

* `clustering.info.tsv` is a DataFrame listing how each gene and scaffold is clustered

#### Intermediate files

* `all.hits` is the parsed HMM hits [maybe garbage?]

* `cluster.txt` is the clustering file produced by VSEARCH

* Each input file will have a file ending in `.faa`, `.genes`, and `.hits`. These are all identified amino-acid sequences, nucleotide sequences, and raw HMM results (respectively)

## Dependencies

* [VSEARCH](https://github.com/torognes/vsearch)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [HMMER](http://hmmer.org/)

## Citations

HMMs were downloaded from https://github.com/MadsAlbertsen/multi-metagenome/blob/master/R.data.generation/essential.hmm,
first published in the following work:

Albertsen, M., Hugenholtz, P., Skarshewski, A., Nielsen, K.L., Tyson, G.W., Nielsen, P.H., 2013. Genome sequences of rare, uncultured bacteria obtained by differential coverage binning of multiple metagenomes. Nat. Biotechnol. 31, 533–538. https://doi.org/10.1038/nbt.2579
