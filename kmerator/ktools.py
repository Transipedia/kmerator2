#!/usr/bin/env python3

"""
tools for Kmerator

- mk-transcripts: build a transcriptome and his jellyfish index for a specie
- mk-genome: build a genome and his jellyfish index for a specie --- COMING SOON
- list-species: list all species available from Ensembl API --- COMING SOON
"""

"""
TODO LIST
- [ ] ADD: mk-genome to create genome and jellyfish of genome
- [ ] ADD: list-species to show all species available from Ensembl API
- [ ] ADD: more species (find if specie is available using ebl-species.py)
- [ ] ADD: mk-transcripts: for current release, use date or size or both to not update if it's seems the same
- [ ] ADD: Handle Ctrl + C
- [ ] ADD: Handle alternative to jellyfish, like kmc and kmtricks
"""


import os
import sys
import argparse
import requests
from bs4 import BeautifulSoup
from urllib.request import urlretrieve
import gzip
# import fileinput      # if
import shutil

import info


__appname__   = "ktools"
__shortdesc__ = "tools for kmerator"
__licence__   = info.LICENCE
__version__   = info.VERSION
__author__    = info.AUTHOR



SPECIES = {
    'human': 'homo_sapiens',
    'mouse': 'mus_musculus',
    'zebrafish': 'danio_rerio'
}



def main():
    """ Function doc """
    args = usage()

    if args.command == 'mk-transcripts':
        MkTranscripts(args)


class MkTranscripts:
    """ Class doc """

    def __init__(self, args):
        """ Class initialiser """

        ### select fasta file to download
        base_url = "http://ftp.ensembl.org/pub/"
        release = f"release-{args.release}/fasta" if args.release else "current_fasta"
        specie = SPECIES[args.specie]

        downloaded_fasta = []
        temp_fasta_files = []

        print(f"{__appname__} {args.command} proceed, please wait...")
        ### download fasta gzipped files
        for item in ('cdna', 'ncrna'):

            ### get url links for cDNA and ncRNA fasta files
            url = os.path.join(base_url, release, specie, item)
            file_name = self.get_link(url, item)
            link = os.path.join(url, file_name)

            ### Download fasta files
            fasta_path = self.wget_fasta(args, link, file_name)
            downloaded_fasta.append(fasta_path)

            ### check fasta file
            self.check_fasta(fasta_path)

            ### create a temp file, with alternate chromosome removed
            temp_fasta_path = self.filtered_fasta(args, fasta_path)
            temp_fasta_files.append(temp_fasta_path)

        ### concatene filtered cDNA and ncRNA (and remove temp fasta files)
        transcriptome_fa = self.mk_transcriptome(temp_fasta_files)

        ### make index of transcriptome (kmc or jellyfish)
        self.mk_index(args, transcriptome_fa)

        ### quit properly
        self.exit_gracefully(args, downloaded_fasta, temp_fasta_files)


    def get_link(self, base_url, item):
        if item == 'cdna':
            pattern = 'cdna.all.fa.gz'
        elif item == 'ncrna':
            pattern = 'ncrna.fa.gz'
        try:
            response = requests.get(base_url, timeout=10)
        except requests.exceptions.ConnectionError:
            sys.exit("Error: unable to connect to Ensembl API.")
        if response.ok:
            response_text = response.text
        else:
            return response.raise_for_status()
        soup = BeautifulSoup(response_text, 'html.parser')
        href = [link.get('href') for link in soup.find_all('a') if pattern in link.get('href')][0]
        return href


    def wget_fasta(self, args, link, file_name):
        release = args.release or 'current'
        gz_file = file_name.split('.')
        gz_file = f"{'.'.join(gz_file[0:2])}.{release}.{'.'.join(gz_file[2:])}"
        # ~ print(link)
        fasta_path = os.path.join(args.output, gz_file)
        ### check if Ensembl fasta files already exists
        if os.path.isfile(fasta_path):
            print(f"Notice: {gz_file!r} already exists, it will no be downloaded.")
            return fasta_path
        os.makedirs(args.output, exist_ok=True)
        urlretrieve(link, fasta_path)
        return fasta_path


    def check_fasta(self, fasta_file):
        try :
            with gzip.open(fasta_file) as fh:
                first_line = fh.readline().rstrip().decode()
                if not first_line.startswith(">"):
                    sys.exit("{}Error: Are you sure {} is a fasta file ?{}".format(bcolors.RED, fasta, bcolors.END))
                return first_line
        except FileNotFoundError:
            sys.exit("{}Error: File '{}' not found.{}".format(bcolors.RED, fasta, bcolors.END))


    def filtered_fasta(self, args, fasta_file_path):
        ## Handle sequences
        sequences = []
        seq = []
        current_header = ""
        end_head =  '\n'            # '\t' if args.tsv else '\n'
        sep_seq  =  '\n'            # '' if args.uniq or args.tsv else '\n'
        with gzip.open(fasta_file_path, 'rt') as fh:
            previous_header = fh.readline()
            for line in fh:
                if line.rstrip():
                    if line.startswith('>'):
                        current_header = line
                        if not 'CHR_' in previous_header.split()[2]:
                            sequences.append(f"{previous_header}{''.join(seq)}\n")
                        seq = []
                        previous_header = current_header
                    else:
                        seq.append(line.rstrip())
        ### last fasta sequence is not printed in the loop
        if not 'CHR_' in previous_header.split()[2]:
            sequences.append(f"{previous_header}{''.join(seq)}\n")
        ### write temp file
        fasta_file_basename = ".".join(os.path.basename(fasta_file_path).split('.')[:-2])
        temp_fasta_path = os.path.join(args.output, f'{fasta_file_basename}-altCHR.fa')
        with open(temp_fasta_path, 'w') as fh:
            for sequence in sequences:
                fh.write(sequence)
        ### check if temparry file done
        if not os.path.isfile(temp_fasta_path):
            sys.exit(f"Error: temporary {temp_fasta_path!r} not found.")
        # ~ print(f"{temp_fasta_path!r} succesfully created")
        return temp_fasta_path


    def output_files_exists(self, args, temp_fasta_files):
        """ Function doc """
        if not os.isdir(args.output):
            return False
        release = args.release or 'current'
        fasta = SPECIE[args.specie]
        transcriptome_fa = f"{'.'.join(temp_fasta_files[0].split('.')[:3])}.cdna+ncrna-altchr.fa"
        if os.path.isfile(transcriptome_fa):
            print(f"Notice: {transcriptome_fa!r} already exists, it will no be downloaded.")
            return True
        return False


    def mk_transcriptome(self, temp_fasta_files):
        ### check if transcriptom already exists
        transcriptome_fa = f"{'.'.join(temp_fasta_files[0].split('.')[:3])}.cdna+ncrna-altchr.fa"
        if os.path.isfile(transcriptome_fa):
            print(f"Notice: {transcriptome_fa!r} already exists, it will be overwritten.")

        ### concatene cDNA and ncRNA
        print(f"creating transcriptome {os.path.basename(transcriptome_fa)!r}")
        with open(transcriptome_fa, 'w') as outfile:
            for fasta in temp_fasta_files:
                with open(fasta, 'r') as infile:
                    # ~ outfile.write(infile.read())        # small files
                    for line in infile:                     # big files
                        outfile.write(line)
        ''' other method
        with open(transcriptome_fa, 'w') as fout, fileinput.input(temp_fasta_files, 'rb') as fin:
            for line in fin:
                fout.write(line)
        '''

        return transcriptome_fa


    def mk_index(self, args, transcriptome_fa):
        """ Function doc """
        ### Select best tool to indexing transcriptome
        tools = ['_kmc', 'jellyfish']
        num = None
        for i,bin in enumerate(tools):
            if shutil.which(bin):
                num = i
                break
        if not num:
            print(f"Warning: no tool found to index the transcriptome")
            return None

        ### Define command line according to the tool used
        ''' in the future (python 3.10)
        match tools[num]:
            case 'kmc':
                print(f"prefered tool: {os.path.basename(tool)}")
            case 'jellyfish':
                cmd = f"{tool} count -m 31 -s 100000 {transcriptome_fa} -o {transcriptome_idx}"
        '''
        tool = tools[num]
        ## kmc case
        if tool == 'kmc': print(f"prefered tool: {os.path.basename(tool)}")
        ## jellyfish case
        elif tool == 'jellyfish':
            transcriptome_idx = f"{os.path.splitext(transcriptome_fa)[0]}.jf"
            cmd = f"{tool} count -t {args.thread} -m 31 -s 100000 {transcriptome_fa} -o {transcriptome_idx}"

        ### Build index
        print(f"Build index with {tool!r}, please wait...")
        ret = os.system(cmd)


    def exit_gracefully(self, args, downloaded_fasta, temp_fasta_files):
        ### remove temp fasta files
        if not args.keep:
            for file in downloaded_fasta:
                os.remove(file)
            for file in temp_fasta_files:
                os.remove(file)



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
                        help = "output directory (default: .)",
                        default=".",
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
                        help="Encode version of transcriptome to dowload (default: current)",
                        )
    parser_mk_tcrptome.add_argument('-t', '--thread',
                        type=int,
                        help="number of trhead (defaul: 1)",
                        default=1,
                       )
    parser_mk_tcrptome.add_argument('-k', '--keep',
                        action="store_true",          # boolean
                        help="Keep intermediates files",
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
