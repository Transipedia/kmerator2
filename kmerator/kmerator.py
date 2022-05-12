#!/usr/bin/env python3

"""
Decomposition of transcript or gene sequences and extraction of their specific k-mers
"""


import sys
import os
import requests
import shutil
import getpass
from datetime import datetime
import signal
from functools import partial

import info
import sequences
from sequences import Sequences
from utils import usage, checkup_args, Color
from kmerize import SpecificKmers
from ensembl import Ensembl



def main():
    ### Handle arguments
    args = usage()
    if args.debug: print(f"{'-'*9}\n{Color.YELLOW}Args: {args}{Color.END}")

    ### Handle control C
    global sigint_handler
    sigint_handler = partial(sigint_handler, args=args)
    signal.signal(signal.SIGINT, sigint_handler)

    ### check options
    checkup_args(args)

    ### some variables
    report = {'aborted': [], 'done': [], 'multiple': []}
    transcriptome_dict = {}
    best_transcripts = {}                       # when genes/transcripts annotated (--selection)
    unannotated_transcripts = []                # when transcripts are unannotated (--fasta-file)

    ### when --selection option is set
    if args.selection:
        ### Load transcriptome as dict (needed to build sequences and to found specific kmers
        print(f" 🧬 Load transcriptome.")
        transcriptome_dict = ebl_fasta2dict(args.transcriptome)
        ### get canonical transcripts using Ensembl API
        print(f" 🧬 Fetch some information from Ensembl API.")
        ensembl = Ensembl(args, report)
        best_transcripts = ensembl.get_ENST()
        ### Build sequence using provided transcriptome
        print(f" 🧬 Build sequences.")
        # ~ build_sequences(args, report, best_transcripts, transcriptome_dict)
        # ~ sequences.build(args, report, best_transcripts, transcriptome_dict)
        seq = Sequences(args, report, best_transcripts, transcriptome_dict)
        seq.build()

    ### when --fasta-file option is set
    else:
        print(f" 🧬 Build sequences.")
        # ~ sequences.build(args, report, unannotated_transcripts)
        seq = Sequences(args, report, unannotated_transcripts)
        seq.build()

    ### get specific kmer with multithreading
    print(f" 🧬 Extract specific kmers, please wait..")
    kmers = SpecificKmers(args, report, transcriptome_dict, best_transcripts, unannotated_transcripts)

    ### Concatene results
    merged_results(args)

    ### show some final info in the prompt
    show_info(report)

    ### set markdown report
    markdown_report(args, report)

    ### ending
    exit_gracefully(args)


def ebl_request(report, item, url, headers):
    try:
        r = requests.get(url, headers=headers)
    except requests.ConnectionError as err:
        sys.exit(f"{Color.RED}\n Error: Ensembl is not accessible or not responding.{Color.END}")
    r = r.json()
    if not r:
        report['aborted'].append(f"{item} not found from Ensembl API.")
        return None
    if 'error' in r:
        report['aborted'].append(f"{r[error]}.")
        return None
    return r


def ebl_fasta2dict(fasta_file):
    '''
    Convertion from Ensembl fasta file to dict
    It keeps only the transcript name of the headers, without version number
    '''
    ### controls
    with open(fasta_file) as fh:
        first_line = fh.readline()
        if not first_line.startswith('>'):
            sys.exit(f"{Color.RED}Error: {os.path.basename(args.fasta_file.name)!r} does not seem to be in fasta format.")
    ### compute file as dict
    fasta_dict = {}
    with open(fasta_file) as fh:
        seq = ""
        old_desc, new_desc = "", ""
        for line in fh:
            if line[0] == ">":
                line = line.split()
                if len(line) > 6:
                    gene_name = line[6].split(':')[1]                    # gene symbol
                else:
                    gene_name = line[3].split(':')[1].split('.')[0]      # ENSG
                transcript_name = line[0].split('.')[0].lstrip('>')
                new_desc = f"{gene_name}:{transcript_name}"
                # ~ new_desc = transcript_name
                if old_desc:
                    fasta_dict[old_desc] = seq
                    seq = ""
                old_desc = new_desc
            else:
                seq += line.rstrip()
    fasta_dict[old_desc] = seq
    return fasta_dict


def merged_results(args):
    if not os.path.isdir(os.path.join(args.output, 'tags')):
        return None
    for item in ['tags', 'contigs']:
        files = os.listdir(os.path.join(args.output, item))
        if files:
            merged_file = os.path.join(args.output, f"{item}-merged.fa")
            with open(merged_file,'wb') as mergefd:
                for file in files:
                    with open(os.path.join(args.output, item, file),'rb') as fd:
                        shutil.copyfileobj(fd, mergefd)


def show_info(report):
    ### show some final info in the prompt
    print(f"{Color.CYAN}\n Done ({len(report['done'])}):")
    for i,mesg in enumerate(report['done']):
        if i == 15:
            print("  - ... (more responses)")
            break
        print(f"  - {mesg}")


    if report['multiple']:
        print(f"{Color.BLUE}\n Multiple responses from Ensembl API ({len(report['multiple'])}):")
        for i,mesg in enumerate(report['multiple']):
            if i == 15:
                print("  - ... (more responses)")
                break
            for k,v in mesg.items():
                print(f"  - {k}: {' '.join(v)}")

    if report['aborted']:
        print(f"{Color.PURPLE}\n Aborted ({len(report['aborted'])}):")
        for i,mesg in enumerate(report['aborted']):
            if i == 15:
                print("  - ... (more responses)")
                break
            print(f"  - {mesg}")

    print(f"{Color.END}")


def markdown_report(args, report):
    with open(os.path.join(args.output, 'report.md'), 'w') as fh:
        fh.write('# kmerator report\n')
        fh.write(f"*date: {datetime.now().strftime('%Y-%m-%d %H:%M')}*  \n")
        fh.write(f'*login: {getpass.getuser()}*\n\n')
        fh.write(f"**kmerator version:** {info.VERSION}\n\n")
        command = ' '.join(sys.argv).replace(' -', ' \\\n  -')
        fh.write(f"**Command:**\n\n```\n{command}\n```\n\n")
        fh.write(f"**Working directory:** `{os.getcwd()}`\n\n")
        fh.write(f"**Jellyfish transcriptome used:** `{args.jellyfish_transcriptome}`\n\n")
        if report['done']:
            fh.write(f"**Genes/transcripts succesfully done ({len(report['done'])})**\n\n")
            for mesg in report['done']:
                fh.write(f"- {mesg}\n")
        if report['multiple']:
            fh.write(f"\n**Multiple Genes returned for one given by Ensembl API ({len(report['multiple'])})**\n\n")
            for mesg in report['multiple']:
                for k,v in mesg.items():
                    fh.write(f"- {k}: {' '.join(v)}\n")
        if report['aborted']:
            fh.write(f"\n\n**Genes/transcript missing ({len(report['aborted'])})**\n\n")
            for mesg in report['aborted']:
                fh.write(f"- {mesg}\n")


def sigint_handler(signal, frame, args):
    exit_gracefully(args)


def exit_gracefully(args):
    if not args.keep:
        shutil.rmtree(os.path.join(args.output,'tags'), ignore_errors=True)
        shutil.rmtree(os.path.join(args.output,'contigs'), ignore_errors=True)
        shutil.rmtree(os.path.join(args.output,'indexes'), ignore_errors=True)
        shutil.rmtree(os.path.join(args.output,'sequences'), ignore_errors=True)
    sys.exit(0)


""" _find_longest_variant
def _find_longest_variant(args, gene_name, transcriptome_dict):
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Finding the longest variant for the gene {gene_name}.{Color.END}")
    variants_dict = { k:len(v) for (k,v) in transcriptome_dict.items() if k.startswith(f"{gene_name}:")}
    # ~ print(*[k for k in variants_dict], sep='\n')
    nb_variants = len(variants_dict)
    if args.verbose: print(f"{Color.YELLOW}Number of variants: {nb_variants}")
    longest_variant = None
    length = 0
    for k,v in variants_dict.items():
        if v > length:
            length = v
            longest_variant = ':'.join(k.split(':')[1:2])
    # ~ print(f"{longest_variant = }")
    return longest_variant
"""


if __name__ == '__main__':
    main()
