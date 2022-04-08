#!/usr/bin/env python3

"""
tools for Kmerator

- mk-transcriptome: build transcriptome for a specie

"""


import os
import sys
import argparse
import requests
from bs4 import BeautifulSoup
from urllib.request import urlretrieve


__appname__   = "k-tools"
__shortdesc__ = "tools for kmerator"
__licence__   = "GPL"
__version__   = "0.1.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"



SPECIES = {
    'human': 'homo_sapiens',
    'mouse': 'mus_musculus',
    'zebrafish': 'danio_rerio'
}



def main():
    """ Function doc """
    args = usage()
    print(args)
    ### select fasta file to download
    base_url = "http://ftp.ensembl.org/pub/"
    release = f"release-{args.release}/fasta" if args.release else "current_fasta"
    specie = SPECIES[args.specie]
    ### download fasta gzipped files
    for item in ('cdna', 'ncrna'):
        ### get url links for cDNA and ncRNA fasta files
        url = os.path.join(base_url, release, specie, item)
        file_name = get_link(url, item)
        link = os.path.join(url, file_name)
        ### Download fasta files
        fasta_file = wget_fasta(args, link, file_name)


def files_exists(args):
    """ Function doc """
    if not os.isdir(args.output):
        return None
    release = args.release or 'current'
    fasta = SPECIE[args.specie]
    return None


def get_link(base_url, item):
    if item == 'cdna':
        pattern = 'cdna.all.fa.gz'
    elif item == 'ncrna':
        pattern = 'ncrna.fa.gz'
    response = requests.get(base_url)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()
    soup = BeautifulSoup(response_text, 'html.parser')
    href = [link.get('href') for link in soup.find_all('a') if pattern in link.get('href')][0]
    return href


def wget_fasta(args, link, file_name):
    release = args.release or 'current'
    gz_file = file_name.split('.')
    gz_file = f"{'.'.join(gz_file[0:2])}.{release}.{'.'.join(gz_file[2:])}"
    print(link)
    dst = os.path.join(args.output, gz_file)
    ### check if final file already exists
    if os.path.isfile(dst):
        print(f"{gz_file!r} already exists, it will no be downloaded.")
        return dst
    os.makedirs(args.output, exist_ok=True)
    urlretrieve(link, dst)
    return dst



def remove_chr_(args):
    pass


def mk_transcriptome(args):
    pass


def jellyfish(args):
    """ Function doc """
    pass



def usage():
    """
    Help function with argument parser.
    https://docs.python.org/3/howto/argparse.html?highlight=argparse
    """
    doc_sep = '=' * min(55, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     prog=__appname__,)
    global_parser = argparse.ArgumentParser(add_help=False)
    subcmd_parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers()
    subparsers.dest = 'command'

    ### OPTION
    global_parser.add_argument('-o', '--output',
                        help = "output directory (default: output)",
                        default="output",
                       )
    parser_mk_tcrptome = subparsers.add_parser('mk-transcripts',
                        parents = [global_parser],
                        help="make transcriptome",
                        )
    parser_mk_tcrptome.add_argument('-s', '--specie',
                        type=str,
                        help="human or mouse (default: human)",
                        default="human",
                        choices=SPECIES,
                        )
    parser_mk_tcrptome.add_argument('-r', '--release',
                        help="Encode version of transcriptome to dowload",
                        )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()
