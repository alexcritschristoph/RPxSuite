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

def main(args):

    hmm_file = os.path.dirname(__file__) + '/helper_files/essential.hmm'
    hmm_name= args.gene

    args.output = args.output.rstrip("/") + "/"
    # if not os.path.isdir(args.output):
    #     os.system("mkdir " + args.output)
    # else:
    #     print("{0} already exists- either remove or choose another output name".format(args.output))
    #     sys.exit()
    #
    # ## MAIN LOOP: CALL PRODIGAL, HMMSEARCH ON EACH FILE
    # for assembly in args.input:
    #
    #     ## call prodigal on contigs
    #     print("Running prodigal on " + assembly)
    #
    #     ## figure out basename of prodigal files
    #     base = args.output + os.path.basename(assembly)
    #
    #     cmd = ["prodigal", '-i', assembly, '-p', "meta", '-d', base + ".genes", '-a', base + ".faa"]
    #     print(' '.join(cmd))
    #     process = subprocess.Popen(cmd, stdout=DEVNULL).wait()
    #     ## call HMMscan on marker genes
    #     print("Running HMMSearch")
    #     if args.score_cutoff == 'cut_tc':
    #         cmd = ["hmmsearch", '--cpu', '6', "--cut_tc", '--tblout', base  + ".hits", hmm_file, base + ".faa"]
    #     elif args.score_cutoff == 'cut_nc':
    #         cmd = ["hmmsearch", '--cpu', '6', "--cut_nc", '--tblout', base  + ".hits", hmm_file, base + ".faa"]
    #     elif args.score_cutoff == 'cut_ga':
    #         cmd = ["hmmsearch", '--cpu', '6', "--cut_ga", '--tblout', base  + ".hits", hmm_file, base + ".faa"]
    #
    #     process = subprocess.Popen(cmd, stdout=DEVNULL).wait()

    ## PART 2
    ## Make FASTA file of hits
    print("Reading hits...")
    f_out = open(args.output + "all.hits", 'w+')
    all_hits = []
    for hit in glob.glob(args.output + "*.hits"):
        f = open(hit)
        hmm = hit.split("/")[-1].split(".")[0]
        for line in f.readlines():
                if not line.startswith("#") and line.split()[2] == hmm_name:
                        hit = line.split()[0]
                        scaf = "_".join(hit.split("_")[:-1])
                        all_hits.append(hit)
    print(all_hits)
    if args.nucleotide:
        suffix = "*.genes"
    else:
        suffix = "*.faa"

    for fn in glob.glob(args.output + suffix):
        for record in SeqIO.parse(fn, "fasta"):
            if record.id in all_hits:
                f_out.write(">" + record.id + "\n")
                f_out.write(str(record.seq) + "\n")
    f_out.close()

    ## Cluster protein hits
    print("Running VSEARCH")
    cmd = ["vsearch", "--cluster_fast", args.output + "all.hits", "--id", str(args.id), "--centroids", args.output + "centroids.fasta", "--uc", args.output + "clusters.txt"]
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= """
        finds marker genes and clusters them and such.\n
        Output: outdir_dir/clusters.txt - vsearch clustering info. outdir/centroids.txt - centroid sequences. outdir/longest.contig.fasta - longest contig sequences for each cluster.\n
        usage: rpX.py *.fasta -g rps3
        """, formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument('input', nargs="+",  help="path to directories, with an asterix. eg. *.fasta, or ./assemblies/*.fasta")
    parser.add_argument("-g", "--gene", action="store", default="Ribosomal_S9", \
        help='Marker gene to look for. Options: PGK,Ribosomal_L23,Ribosomal_L5,Ribosomal_L3,Ribosomal_L6,Ribosomal_S17,Ribosomal_S9,Ribosomal_S8,Ribosomal_S11,Ribosomal_S13,Ribosomal_L10,Ribosomal_L4,tRNA-synt_1d,GrpE,Methyltransf_5,TIGR00001,TIGR00002,TIGR00009,TIGR00012,TIGR00019,TIGR00029,TIGR00043,TIGR00059,TIGR00060,TIGR00061,TIGR00062,TIGR00064,TIGR00082,TIGR00086,TIGR00092,TIGR00115,TIGR00116,TIGR00152,TIGR00158,TIGR00165,TIGR00166,TIGR00168,TIGR00234,TIGR00337,TIGR00344,TIGR00362,TIGR00388,TIGR00389,TIGR00392,TIGR00396,TIGR00408,TIGR00409,TIGR00414,TIGR00418,TIGR00420,TIGR00422,TIGR00435,TIGR00436,TIGR00442,TIGR00459,TIGR00460,TIGR00468,TIGR00471,TIGR00472,TIGR00487,TIGR00496,TIGR00575,TIGR00631,TIGR00663,TIGR00775,TIGR00810,TIGR00855,TIGR00922,TIGR00952,TIGR00959,TIGR00963,TIGR00964,TIGR00967,TIGR00981,TIGR01009,TIGR01011,TIGR01017,TIGR01021,TIGR01024,TIGR01029,TIGR01030,TIGR01031,TIGR01032,TIGR01044,TIGR01049,TIGR01050,TIGR01059,TIGR01063,TIGR01066,TIGR01067,TIGR01071,TIGR01079,TIGR01164,TIGR01169,TIGR01171,TIGR01391,TIGR01393,TIGR01632,TIGR01953,TIGR02012,TIGR02013,TIGR02027,TIGR02191,TIGR02350,TIGR02386,TIGR02387,TIGR02397,TIGR02432,TIGR02729,TIGR03263,TIGR03594')
    parser.add_argument("-i", "--id", action="store", default=0.97, \
        help='Marker gene percent identity cutoff to use when clustering. ')
    parser.add_argument("-o", "--output", action="store", default="./marker_genes/", \
        help='Output directory')
    parser.add_argument("-s", "--score_cutoff", action="store", default="cut_ga",  \
        help='An HMM score threshold to use - cut_ga, cut_nc, or cut_tc.')
    parser.add_argument("--nucleotide", dest='nucleotide', action='store_true', \
        help='Choose to cluster by nucleotide sequence instead of protein sequence.')

    parser.set_defaults(nucleotide=False)
    args = parser.parse_args()
    main(args)
