# RPxSuite
RPxSuite is a dedicated, stand-alone, tool for the identification and clustering of marker genes in shotgun metagenomic assemblies.

The publication supporting this tool is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/647511v1).

Briefly, this tool aims to automate a portion of what was done in the methods section of Diamond et al. 2019:
https://www.nature.com/articles/s41564-019-0449-y#Sec12

```
rpS3 marker sequences were identified across all metagenomes using a custom hidden Markov model (HMM) based on an alignment of rpS3 sequences from the tree of life data set from ref. 51. Briefly, all rpS3 sequences provided in ref. 51 were initially filtered to remove Eukaryotic sequences. Sequences were then clustered at 90% ID using USEARCH with the following parameters: usearch -cluster_fast rpS3_sequences.faa -sort length -id 0.90 -maxrejects 0 -maxaccepts 0 -centroids rpS3_sequences_NR90.faa. The non-redundant sequences were then filtered to remove sequences <200 amino acids in length with pullseq. The resulting 2,249 sequences were aligned using muscle52 and an HMM was constructed from the alignment using HMMER3 with default parameters53. The HMM was benchmarked against the Uniprot reference proteomes database, and it was determined that rpS3 sequences could be confidently identified above a cutoff HMM alignment score of 40.

Across all metagenomes we identified a total of 10,159 rpS3 sequences that passed our HMM score threshold of 40. We clustered these sequences at 99% ID using USEARCH to obtain groups that roughly equate to species. We refer to these as species groups (SGs). The following USEARCH options were used: -cluster_fast all_rpS3.fa -sort length -id 0.99 -maxrejects 0 -maxaccepts 0 -centroids all_rpS3_centroids.faa. Subsequently we identified the longest contig in each rpS3 protein cluster to serve as a mapping target for abundance quantification of each SG (Supplementary Table 3 and Supplementary Data 2).
```

## Installation with pip
```
$ git clone https://github.com/alexcritschristoph/RPxSuite.git

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
