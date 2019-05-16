# RPxSuite
RPxSuite is a dedicated, stand-alone, tool for the identification and clustering of marker genes in shotgun metagenomic assemblies.

The publication supporting this tool is available on [bioRxiv](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

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

## Dependencies

* [VSEARCH](https://github.com/torognes/vsearch)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [HMMER](http://hmmer.org/)

## Citations

HMMs were downloaded from https://github.com/MadsAlbertsen/multi-metagenome/blob/master/R.data.generation/essential.hmm,
first published in the following work:

Albertsen, M., Hugenholtz, P., Skarshewski, A., Nielsen, K.L., Tyson, G.W., Nielsen, P.H., 2013. Genome sequences of rare, uncultured bacteria obtained by differential coverage binning of multiple metagenomes. Nat. Biotechnol. 31, 533–538. https://doi.org/10.1038/nbt.2579
