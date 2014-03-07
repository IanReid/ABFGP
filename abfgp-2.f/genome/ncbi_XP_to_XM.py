__author__ = 'ian'

import sys, os

this_dir = os.path.dirname(__file__)
src = os.path.dirname(this_dir)
sys.path.append(src)
import argparse
from Bio import Entrez

DESCRIPTION = 'Find the XM_ transcript names corresponding to XP_ protein names'
VERSION = '0.1'
BATCH_SIZE = 50  # Number of names to process in one request


def get_args():
    argparser = argparse.ArgumentParser(description=DESCRIPTION)
    # standard options
    argparser.add_argument('--version', action='version', version='%(prog)s' + VERSION)
    argparser.add_argument('--verbose', '-v', action='count', default=0,
                           help='Omit to see only fatal error messages; -v to see warnings; -vv to see warnings and progress messages')
    # options to customize
    argparser.add_argument('--in', '-i', dest='input', type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                           help='Path to the input file; if omitted or -, input is read from stdin')
    argparser.add_argument('--out', '-o', type=argparse.FileType('w'), nargs='?', default=sys.stdout,
                           help='Path to the output file; if omitted or -, output is written to stdout')
    return argparser.parse_args()


def ncbi_XP_to_XM(batch):
    """
    :param: batch
    :type: String of protein IDs joined with commas
    :return:
    :rtype: String of mRNA IDs joined with newlines
    """
    search_result = Entrez.read(Entrez.esearch(db='protein',term=batch,usehistory='y'))
    if search_result['Count'] == '0':
        raise ValueError('No protein database matches for ' + batch)
    qkey =int(search_result['QueryKey'])
    wenv = search_result['WebEnv']
    link_result = Entrez.read(Entrez.elink(db='nuccore',dbfrom='protein',WebEnv=wenv,query_key=qkey,linkname='protein_nuccore_mrna',cmd='neighbor_history'))[0]
    qkey = int(link_result['LinkSetDbHistory'][0]['QueryKey'])
    wenv = link_result['WebEnv']
    final_result = Entrez.efetch(db='nuccore',query_key=qkey,WebEnv=wenv,rettype='acc',retmode='text').read()
    return final_result


if __name__ == '__main__':
    args = get_args()
    start = 0
    p2t = {}
    protein_ids = []
    for line in args.input:
        protein_ids.extend(line.strip().split(','))

    while start < len(protein_ids):
        slice = protein_ids[start:start+BATCH_SIZE]
        batch = ','.join(slice)
        transcript_batch = ncbi_XP_to_XM(batch)
        transcript_ids = transcript_batch.splitlines()
        for p, t in zip(slice,transcript_ids):
            p2t[p] = t
            print >> args.out, '%s\t%s' % (p, t)
        start += BATCH_SIZE

    args.out.close()
    print >> sys.stderr, sys.argv[0], 'done.'
