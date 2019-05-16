#!/usr/bin/env python


import os
import sys
import glob
import argparse
import subprocess
from Bio import SeqIO
from subprocess import DEVNULL
from collections import defaultdict

def main(args):

    hmm_file = os.path.dirname(__file__) + '/helper_files/essential.hmm'
    hmm_name= args.gene

    args.output = args.output.rstrip("/") + "/"
    if not os.path.isdir(args.output):
        os.system("mkdir " + args.output)
    else:
        print("{0} already exists- either remove or choose another output name".format(args.output))
        sys.exit()

    ## MAIN LOOP: CALL PRODIGAL, HMMSEARCH ON EACH FILE
    for assembly in args.input:

        ## call prodigal on contigs
        print("Running prodigal on " + assembly)

        cmd = ["prodigal", '-i', assembly, '-p', "meta", '-d', args.output + assembly + ".genes", '-a', args.output + assembly + ".faa"]
        process = subprocess.Popen(cmd, stdout=DEVNULL).wait()
        ## call HMMscan on marker genes
        print("Running HMMSearch")
        if args.score_cutoff == 'cut_tc':
            cmd = ["hmmsearch", '--cpu', '6', "--cut_tc", '--tblout', args.output + assembly  + ".hits", hmm_file, args.output + assembly + ".faa"]
        elif args.score_cutoff == 'cut_nc':
            cmd = ["hmmsearch", '--cpu', '6', "--cut_nc", '--tblout', args.output + assembly  + ".hits", hmm_file, args.output + assembly + ".faa"]
        elif args.score_cutoff == 'cut_ga':
            cmd = ["hmmsearch", '--cpu', '6', "--cut_ga", '--tblout', args.output + assembly  + ".hits", hmm_file, args.output + assembly + ".faa"]

        process = subprocess.Popen(cmd, stdout=DEVNULL).wait()

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

    ## Get centroid information
    centroids = {}
    for record in SeqIO.parse(args.output + "centroids.fasta", "fasta"):
        centroids["_".join(record.id.split("_")[:-1])] = record.id

    f = open(args.output + "clusters.txt")
    clusters = defaultdict(list)
    for line in f.readlines():
        cluster = line.split()[1]
        prot = line.split()[-2]
        scaf = "_".join(prot.split("_")[:-1])
        clusters[cluster].append(scaf)
    f.close()

    ## Get scaffold lengths
    s2l = {}
    seqs = {}
    for fn in args.input:
        for record in SeqIO.parse(fn, "fasta"):
            s2l[record.id] = len(record.seq)
            seqs[record.id] = record.seq

    ## Get largest contig
    f_long = open(args.output + "longest.contigs.fasta", "w+")
    print("Cluster Number\tLongest Scaffold\tCentroid Scaffold\tGenesInScaffold")
    for cluster in clusters:
        largest = 0
        largest_name = 'NA'
        centroid = 'NA'
        for scaf in clusters[cluster]:
            if s2l[scaf] > largest:
                largest = s2l[scaf]
                largest_name = scaf
            if scaf in centroids:
                centroid = centroids[scaf]
        print(cluster + "\t" + largest_name + "\t" + centroid + "\t" + str(len(clusters[cluster])))
        f_long.write(">" + cluster + "_" + largest_name + ":" + centroid + "." + str(len(clusters[cluster])) + "\n")
        f_long.write(str(seqs[largest_name]) + "\n")

    f_long.close()

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
