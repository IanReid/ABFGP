__author__ = 'ian'

import sys, os

this_dir = os.path.dirname(__file__)
src = os.path.dirname(this_dir)
sys.path.append(src)
import argparse

DESCRIPTION = 'What this program does'
VERSION = '0.1'


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


if __name__ == '__main__':
    args = get_args()

    args.out.close()
    print >> sys.stderr, sys.argv[0], 'done.'
