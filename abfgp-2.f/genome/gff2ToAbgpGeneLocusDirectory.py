__author__ = 'ian'

import sys, os

this_dir = os.path.dirname(__file__)
import argparse
from Bio import SeqIO
import re
from itertools import groupby

name_pat = re.compile('name "([^"]+)";')

DESCRIPTION = 'Create a set of ABFGP GeneLocusDirectories for the protein coding genes annotated in a GFF2 file'
VERSION = '0.1'
FLANK_LENGTH = 1000

def get_args():
    argparser = argparse.ArgumentParser(description=DESCRIPTION)
    # standard options
    argparser.add_argument('--version', action='version', version='%(prog)s' + VERSION)
    argparser.add_argument('--verbose', '-v', action='count', default=0,
                           help='Omit to see only fatal error messages; -v to see warnings; -vv to see warnings and progress messages')
    # options to customize
    argparser.add_argument('--in', '-i', dest='input', type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                           help='Path to the input GFF3 annotation file; if omitted or -, input is read from stdin')
    argparser.add_argument('--genome','-g',required=True,help='Genome sequence file in FASTA format')
    argparser.add_argument('--outdir','-o',required=True,help='Path to output root directory')
    argparser.add_argument('--tag','-t',dest='organism_tag',required=True,help='Unique identifier for this genome')
    return argparser.parse_args()

def get_name(line):
    name = None
    match = name_pat.search(line)
    if match:
        name = match.group(1)
    return name

if __name__ == '__main__':
    args = get_args()
    outdir = os.path.join(args.outdir,args.organism_tag)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    genome = SeqIO.to_dict(SeqIO.parse(open(args.genome),'fasta'))

    for genelocusid,gene_lines_iter in groupby(args.input.xreadlines(),key=get_name):
        if genelocusid:
            genelocusdir = os.path.join(outdir,genelocusid)
            if not os.path.isdir(genelocusdir):
                os.makedirs(genelocusdir)
            gff_tuples = [line.strip().split('\t') for line in gene_lines_iter]
            chrom_id = gff_tuples[0][0]
            gene_start = int(gff_tuples[0][3])
            gene_end = int(gff_tuples[0][4])
            for gff_tuple in gff_tuples[1:]:
                gene_start = min(gene_start, int(gff_tuple[3]))
                gene_end = max(gene_end,int(gff_tuple[4]) )
            strand = gff_tuples[0][6]

            # Make dna.fa file
            chrom = genome[chrom_id]
            start = max(1, gene_start - FLANK_LENGTH)
            end = min(len(chrom)-1,gene_end + FLANK_LENGTH)
            sequence = chrom[start-1:end]
            seq_id = '%s.%s' % (args.organism_tag, genelocusid)
            sequence.id = seq_id
            sequence.description = 'Padded genome sequence'
            if strand == '-':
                sequence.seq = sequence.seq.reverse_complement()
                sequence.description += ' reverse-complemented'
            with open(os.path.join(genelocusdir,'%s.dna.fa' % genelocusid),'w') as seqout:
                SeqIO.write(sequence,seqout,'fasta')

            # Make locus.gff file
            with open(os.path.join(genelocusdir,'%s.locus.gff' % genelocusid),'w') as locusout:
                locus_str = '%s\tabgp_locus\tabgp_locus\t%s\t%s\t.\t%s\t.\tgene_id "%s"; abgp_locus "%s"' % \
                            (chrom_id,start, end, strand,seq_id,seq_id)
                locusout.write(locus_str)

            # Make gene.gff file
            with open(os.path.join(genelocusdir,'%s.gene.gff' % genelocusid),'w') as geneout:
                for gff_tuple in gff_tuples:
                    gff_tuple[8] = 'gene_id "%s"; %s' % (genelocusid, gff_tuple[8])
                    print >> geneout, '\t'.join(gff_tuple)

    print >> sys.stderr, sys.argv[0], 'done.'
