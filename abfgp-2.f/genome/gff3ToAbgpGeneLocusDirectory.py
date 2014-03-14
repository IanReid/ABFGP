__author__ = 'ian'

import sys, os

this_dir = os.path.dirname(__file__)
# src = os.path.dirname(this_dir)
# sys.path.append(src)
sys.path.append('/home/ian/python/')
import argparse
from Bio import SeqIO
from transcriptomics.src.GeneModels.gff3Iterator import GFF3Iterator

DESCRIPTION = 'Create a set of ABFGP GeneLocusDirectories for the protein coding genes annotated in a GFF3 file'
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


if __name__ == '__main__':
    args = get_args()
    outdir = os.path.join(args.outdir,args.organism_tag)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    genome = SeqIO.to_dict(SeqIO.parse(open(args.genome),'fasta'))

    for gene in GFF3Iterator(args.input).genes():
        for transcript in gene.get_transcripts():
            genelocusid = transcript.get_ID()
            genelocusdir = os.path.join(outdir,genelocusid)
            if not os.path.isdir(genelocusdir):
                os.makedirs(genelocusdir)

            # Make dna.fa file
            chrom = genome[transcript.get_seqID()]
            start = max(1, transcript.get_start() - FLANK_LENGTH)
            end = min(len(chrom)-1,transcript.get_end() + FLANK_LENGTH)
            sequence = chrom[start-1:end]
            seq_id = '%s.%s' % (args.organism_tag, genelocusid)
            sequence.id = seq_id
            sequence.description = 'Padded genome sequence'
            if transcript.get_strand() == '-':
                sequence.seq = sequence.seq.reverse_complement()
                sequence.description += ' reverse-complemented'
            with open(os.path.join(genelocusdir,'%s.dna.fa' % genelocusid),'w') as seqout:
                SeqIO.write(sequence,seqout,'fasta')

            # Make locus.gff file
            with open(os.path.join(genelocusdir,'%s.locus.gff' % genelocusid),'w') as locusout:
                locus_str = '%s\tabgp_locus\tabgp_locus\t%s\t%s\t.\t%s\t.\tgene_id "%s"; abgp_locus "%s"' % \
                            (transcript.get_seqID(),start, end, transcript.get_strand(),seq_id,seq_id)
                locusout.write(locus_str)

            # Make gene.gff file
            with open(os.path.join(genelocusdir,'%s.gene.gff' % genelocusid),'w') as geneout:
                attributes = '\tgene_id "%s"; name "%s"; transcriptId "%s"' %(seq_id,seq_id,seq_id)
                start_codon = str(transcript.make_start_codon()).rsplit('\t',1)[0] + attributes
                stop_codon = str(transcript.make_stop_codon()).rsplit('\t',1)[0] + attributes

                # ABFGP wants to see CDS lines
                # We assume that if any CDS exon is present, they all are
                # Otherwise we make CDS lines from the exon lines
                CDSpresent = any([exon.get_type() == 'CDS' for exon in transcript.get_exons()])
                CDSstart = CDSend = -1
                if transcript.get_strand() == '-':
                    print >> geneout, stop_codon
                    CDSstart = transcript.get_CDS_stop()
                    CDSend = transcript.get_CDS_start()
                else:
                    print >> geneout,start_codon
                    CDSstart = transcript.get_CDS_stop()
                    CDSend = transcript.get_CDS_start()

                for exon in sorted(transcript.get_exons()):
                    exon_str = str(exon).rsplit('\t',1)[0] + attributes
                    print >> geneout, exon_str
                    if not CDSpresent and exon.get_start() < CDSend and exon.get_end() > CDSstart:
                        cds = exon.clone()
                        cds.set_type('CDS')
                        cds.set_start(max(CDSstart, cds.get_start()))
                        cds.set_end(min(CDSend,cds.get_end()))
                        cds_str = str(cds).rsplit('\t',1)[0] + attributes
                        print >> geneout, cds_str

                if transcript.get_strand() == '-':
                    print >> geneout, start_codon
                else:
                    print >> geneout,stop_codon






    print >> sys.stderr, sys.argv[0], 'done.'
