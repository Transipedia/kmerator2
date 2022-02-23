#!/usr/bin/env python3

'''
les arguments sont une liste de mots
'''

import sys
import argparse

def usage():
    genes = '' if sys.stdin.isatty() else sys.stdin.read().rstrip().split()
    parser = argparse.ArgumentParser()
    parser.add_argument("genes",
                        help="gene list or stdin (example: RUNX1 ENSG00000157764",
                        nargs='*' if genes else '+',
                        default=genes,
                        metavar=('gene'),
                       )
    if len(sys.argv) == 1 and sys.stdin.isatty():
        sys.exit(parser.print_help())
    return parser.parse_args()

args = usage()
print(args)
