#!/usr/bin/env python

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
from subprocess import DEVNULL
from collections import defaultdict

import rpXsuite

def main(args):

    hmm_file = os.path.dirname(rpXsuite.__file__) + '/helper_files/essential.hmm'
    hmm_name= args.gene

    args.output = args.output.rstrip("/") + "/"
    if not os.path.isdir(args.output):
        os.system("mkdir " + args.output)
    else:
        print("{0} already exists- either remove or choose another output name".format(args.output))
        sys.exit()

    ## MAIN LOOP: CALL PRODIGAL, HMMSEARCH ON EACH FILE
    for assembly in args.input:

        ## figure out basename of prodigal files
        base = args.output + os.path.basename(assembly)

        ## call prodigal on contigs
        if not args.prodigal:
            print("Running prodigal on " + assembly)


            cmd = ["prodigal", '-i', assembly, '-p', "meta", '-d', base + ".genes", '-a', base + ".faa"]
            print(' '.join(cmd))
            process = subprocess.Popen(cmd, stdout=DEVNULL).wait()
            prodigal_file_name = base + '.faa'
        else:
            if args.amino_acid:
                prodigal_file_name = args.prodigal
            else:
                ## need to translate from .FNA to .FAA for HMMSEARCH
                print("Translating genes")
                try:
                    f_out = open(base + "_translated_genes.faa", "w+")
                    for record in SeqIO.parse(args.prodigal, 'fasta'):
                        f_out.write(">" + str(record.id) + "\n")
                        f_out.write(str(record.seq.translate()) + "\n")
                    f_out.close()
                except:
                    print("ERROR: Error in translating your prodigal genes. You probably passed a prodigal Protein file when you did not select --amino_acid")
                    sys.exit(1)
                prodigal_file_name = base + "_translated_genes.faa"
        ## call HMMscan on marker genes
        print("Running HMMSearch")
        if args.score_cutoff == 'cut_tc':
            cmd = ["hmmsearch", '--cpu', '6', "--cut_tc", '--tblout', base  + ".hits", hmm_file, prodigal_file_name]
        elif args.score_cutoff == 'cut_nc':
            cmd = ["hmmsearch", '--cpu', '6', "--cut_nc", '--tblout', base  + ".hits", hmm_file, prodigal_file_name]
        elif args.score_cutoff == 'cut_ga':
            cmd = ["hmmsearch", '--cpu', '6', "--cut_ga", '--tblout', base  + ".hits", hmm_file, prodigal_file_name]

        process = subprocess.Popen(cmd, stdout=DEVNULL).wait()

    ## PART 2
    ## Make FASTA file of hits
    print("Reading hits...")
    f_out = open(args.output + "all.hits", 'w+')
    all_hits = []
    for hit in glob.glob(args.output + "*.hits"):
        f = open(hit)
        hmm = hit.split("/")[-1].split(".hits")[0]
        for line in f.readlines():
                if not line.startswith("#") and line.split()[2] == hmm_name:
                        hit = line.split()[0]
                        scaf = "_".join(hit.split("_")[:-1])
                        all_hits.append(hit)
    print(all_hits)
    if args.amino_acid:
        suffix = "*.faa"
    else:
        suffix = "*.genes"

    if not args.prodigal:
        for fn in glob.glob(args.output + suffix):
            for record in SeqIO.parse(fn, "fasta"):
                if record.id in all_hits:
                    f_out.write(">" + record.id + "\n")
                    f_out.write(str(record.seq) + "\n")
    if args.prodigal:
        for record in SeqIO.parse(args.prodigal, "fasta"):
                if record.id in all_hits:
                    f_out.write(">" + record.id + "\n")
                    f_out.write(str(record.seq) + "\n")
    f_out.close()


    ## Cluster protein hits
    print("Running VSEARCH")
    cmd = ["vsearch", "--cluster_fast", args.output + "all.hits", "--id", str(args.id), "--centroids", args.output + "centroids.fasta", "--maxrejects", "0", "--threads", "6", "--uc", args.output + "clusters.txt"]
    print(' '.join(cmd))
    process = subprocess.Popen(cmd, stdout=DEVNULL).wait()

    ## Return cluster information
    print("Reading clustering results")

    Rdb = parse_usearch_clustering(args.output + "clusters.txt")
    Rdb['scaffold'] = ["_".join(prot.split("_")[:-1]) for prot in Rdb['sequence']]
    Rdb['centroid_scaffold'] = ["_".join(prot.split("_")[:-1]) for prot in Rdb['centroid']]

    ## Get scaffold lengths
    s2l = {}
    seqs = {}
    for fn in args.input:
        for record in SeqIO.parse(fn, "fasta"):
            s2l[record.id] = len(record.seq)
            seqs[record.id] = record.seq
    Rdb['length'] = Rdb['scaffold'].map(s2l)

    ## Get largest contig
    f_long = open(args.output + "longest.contigs.fasta", "w+")
    for cluster, db in Rdb.groupby('cluster'):
        cluster = str(cluster)
        leng = str(db.sort_values('length', ascending=False)['length'].tolist()[0])
        largest_name = str(db.sort_values('length', ascending=False)['scaffold'].tolist()[0])
        centroid = str(db.sort_values('length', ascending=False)['cluster'].tolist()[0])
        f_long.write(">" + cluster + "_" + largest_name + ":" + centroid + "." + str(leng) + "\n")
        f_long.write(str(seqs[largest_name]) + "\n")
    f_long.close()

    ## Print info
    Rdb = Rdb.rename(columns={'sequence':'gene', 'centroid':'centroid_gene', 'length':'scaffold_length'})
    Rdb.to_csv(args.output + "clustering.info.tsv", sep='\t', index=False)

def parse_usearch_clustering(loc):
    '''
    From the location of a .uc usearch file, return something like Cdb

    https://www.drive5.com/usearch/manual/cmd_calc_distmx.html
    https://www.drive5.com/usearch/manual/opt_uc.html
    '''
    dtypes = {0:'category', 1:'category', 2:np.int32, 8:'object'}
    ucols = [0,1,2,8]
    Rdb = pd.read_csv(loc, header=None, usecols=ucols,\
            dtype=dtypes, sep='\t')
    table = defaultdict(list)

    # Find the centroids
    sdb  = Rdb[Rdb[0] == 'S']
    shdb = Rdb[Rdb[0].isin(['H', 'S'])]
    for centroid, cdb in sdb.groupby(1):
        cent = cdb[8].tolist()[0].split()[0]
        db = shdb[shdb[1] == centroid]

        #assert len(db) == int(Rdb[2][(Rdb[0] == 'C') & (Rdb[1] == centroid)].tolist()[0])

        for seq in db[8].tolist():
            table['cluster'].append(int(centroid))
            table['members'].append(len(db))
            table['sequence'].append(seq.split()[0])
            table['centroid'].append(cent)

    return pd.DataFrame(table)

def load_genes(describe = False):
    '''
    Load gene thresholds and return g2t (gene -> threshold)

    If describe, print a description of the genes available
    '''
    # Get the thresholds and recoverability
    Rdb = pd.read_csv(os.path.dirname(rpXsuite.__file__) + '/helper_files/SupplementalTable_S4.2.tsv', sep='\t')
    g2t = Rdb.set_index('gene')['overall_threshold'].to_dict()

    if describe:
        # Get the descriptions
        Ddb = pd.read_csv(os.path.dirname(rpXsuite.__file__) + '/helper_files/SupplementalTable_S2.2.tsv', sep='\t')
        g2d = Ddb.sort_values('description', ascending=False).set_index('gene')['description'].to_dict()

        # Add
        Rdb['description'] = Rdb['gene'].map(g2d)

        # Alter
        Rdb = Rdb[~(Rdb['gene'] == '16S')]
        Rdb = Rdb.rename(columns={'overall_threshold':'species_threshold'})
        Rdb = Rdb[['gene', 'species_threshold', 'description']]
        Rdb = Rdb[~Rdb['description'].isna()]

        # Print
        print(Rdb.to_string(index=False, max_rows=1000, line_width=200))

    return(g2t)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= """
        finds marker genes and clusters them and such.\n
        Output: outdir_dir/clusters.txt - vsearch clustering info. outdir/centroids.txt - centroid sequences. outdir/longest.contig.fasta - longest contig sequences for each cluster.\n
        usage: rpX.py *.fasta -g rps3
        """, formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument('input', nargs="+",  help="path to directories, with an asterix. eg. *.fasta, or ./assemblies/*.fasta")
    parser.add_argument("-g", "--gene", action="store", default="Ribosomal_L6", \
        help='Marker gene to look for. Use the argument "--describe_genes" for a list of options')
    parser.add_argument("-i", "--id", action="store", default=0, \
        help='Marker gene percent identity cutoff to use when clustering. By default use the optimal threshold for species-level clustering')
    parser.add_argument("-o", "--output", action="store", default="./marker_genes/", \
        help='Output directory')
    parser.add_argument("-s", "--score_cutoff", action="store", default="cut_ga",  \
        help='An HMM score threshold to use - cut_ga, cut_nc, or cut_tc.')
    parser.add_argument("--amino_acid", dest='amino_acid', action='store_true', \
        help='Choose to cluster by amino_acid sequence instead of nucleotide sequence.')
    parser.add_argument("-p", "--prodigal", action="store", default=None, \
        help='A prodigal predicted proteins file (output by prodigal -d or prodigal -a) - will skip running prodigal if provided. Make sure this file matches your --amino_acid parameter.')
    parser.add_argument("--describe_genes", action='store_true', \
        help='Print the gene options and exit')

    parser.set_defaults(amino_acid=False)
    parser.set_defaults(describe_genes=False)
    args = parser.parse_args()

    # describe genes
    if args.describe_genes == True:
        load_genes(describe=True)
        sys.exit()

    # Get the threshold
    if args.id == 0:
        g2t = load_genes()
        thresh = g2t[args.gene]
        args.id = thresh
        print("ID threshold of {0} being used for gene {1}".format(args.id, args.gene))

    main(args)
