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

import info
from utils import usage, checkup_args, Color
from fasta import fasta2dict
from kmerize import SpecificKmers
from ensembl import Ensembl



BASE_URL = "https://rest.ensembl.org"


def main():
    ### Handle arguments
    args = usage()
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Args:\n{args}{Color.END}")

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
        ensembl = Ensembl(args, report, transcripts)
        best_transcripts = ensembl.get_ENST()
        # ~ best_transcripts = get_ensembl_transcripts(args, report)
        # ~ print(best_transcripts)
        ### Build sequence using provided transcriptome
        print(f" 🧬 Build sequences.")
        build_sequences(args, report, best_transcripts, transcriptome_dict)
    ### when --fasta-file option is set
    else:
        print(f" 🧬 Build sequences.")
        build_sequences(args, report, unannotated_transcripts)

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
    gracefully_exit(args)



def build_sequences(args, report, transcripts, transcriptome_dict=None):
    ''''
    create files for each transcript
    Be careful:
        transcripts == best_transcripts  when genes/transcripts are known
        transcripts == unannotated_transcripts  for unannotated sequences
    '''
    output_seq_dir = os.path.join(args.output, 'sequences')
    removed_transcripts = []
    ### Whith --selection option
    if args.selection:
        ### Abort if no transcripts found
        if not transcripts:
            sys.exit(f"{Color.RED}Error: no sequence found for {args.selection}")
        ### create output directory structure
        os.makedirs(output_seq_dir, exist_ok=True)
        ### Get the sequences and create files for each of them
        for transcript,values in transcripts.items():
            desc = f"{values['symbol']}:{transcript}"
            if desc in transcriptome_dict:
                seq = transcriptome_dict[desc]
                if len(seq) < args.kmer_length:
                    report['warming'].append(f"{desc!r} sequence length < {args.kmer_length} => ignored")
                    continue
                ### create fasta files
                outfile = f"{values['symbol'].replace('.','_')}.{transcript}.fa"[:255].replace(' ', '_').replace('/', '@SLASH@')
                outfile = f"{args.output}/sequences/{outfile}"
                with open(outfile, 'w') as fh:
                    fh.write(f">{values['symbol']}:{transcript}\n{seq}")
            ### When transcript is not found
            else:
                report['aborted'].append(f"{transcript} not found in provided transcriptome (gene: {values['symbol']})")
                removed_transcripts.append(transcript)
            '''
            ### As alternative, fetch sequences with Ensembl API
            ext_ebl = f'/sequence/id/{transcript}?type=cdna;species={args.specie}'
            r = requests.get(BASE_URL+ext_ebl, headers={ "Content-Type": "text/plain"})
            seq = r.text
            '''
        for tr in removed_transcripts:
            transcripts.pop(tr)
    ### Whith --fasta-file option
    else:
        ### read fasta file
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nBuild sequences without transcriptome.\n{Color.END}")
        fastafile_dict = fasta2dict(args.fasta_file)
        ### Abort if dict empy
        if not fastafile_dict:
            sys.exit(f"{Color.RED}Error: no sequence found for {args.fasta_file}")
        ### create output directory structure
        os.makedirs(output_seq_dir, exist_ok=True)
        for desc,seq in fastafile_dict.items():
            outfile = f"{desc.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            if len(seq) < args.kmer_length:
                report['aborted'].append(f"{desc!r} sequence length < {args.kmer_length} => ignored")
                continue
            transcripts.append(desc)
            with open(os.path.join(output_seq_dir, outfile), 'w') as fh:
                fh.write(f">{desc[:79]}\n{seq}")


def get_ensembl_transcripts(args, report):
    ''''
    Works with --selection option,
    - get canonical transcript and symbol name if ENSG is provided
    - get symbol name if ENST is provided return dict as format {ENST: SYMBOL}
    - get canonical transcript if symbol name is provided
    return dict as format {ENST: SYMBOL}
    '''
    ### Define genes/transcripts provided when they are in a file
    if len(args.selection) == 1 and os.path.isfile(args.selection[0]):
        with open(args.selection[0]) as fh:
            args.selection = fh.read().split()
    ### get transcript from Ensembl API
    transcripts = {}                                # the dict to return
    ebl_motifs = ['ENSG', 'ENST']                   # Accepted Ensembl motifs
    ext_symbol = f"/xrefs/symbol/{args.specie}/"    # extension to get ENSG with SYMBOL
    ext_ebl = "/lookup/id/"                         # extension to get info with ENSG/ENST
    headers={ "Content-Type" : "application/json"}  # header for the query
    for item in args.selection:
        ### When ENST is provided, get symbol
        if item.startswith('ENST'):
            url = BASE_URL+ext_ebl+item+"?"
            r = ebl_request(report, item, url, headers=headers)
            if not r: continue
            if  not 'display_name' in r:
                print(f"display name of {item!r} not found from Ensembl API.")
                continue
            transcript = item.split('.')[0]
            symbol = r['display_name'].split('-')[0]
            transcripts[transcript] = {'symbol':symbol, 'level': 'transcript', 'given': item}
        ### When ENSEMBL GENE NAME is provided, get canonical transcript
        elif item.startswith('ENS'):
            url = BASE_URL+ext_ebl+item+"?"
            r = ebl_request(report, item, url, headers=headers)
            if not r: continue
            transcript = r['canonical_transcript'].split('.')[0]
            if 'display_name' in r:
                symbol = r['display_name']
            else:
                symbol = r['id']
            transcripts[transcript] = {'symbol':symbol, 'level': 'gene', 'given': item}
        ### In other cases, item is considered as NAME_SYMBOL
        else:
            url = BASE_URL+ext_symbol+item+"?"
            r = ebl_request(report, item, url, headers=headers)
            if not r:
                continue
            candidates_symbol = []
            for a in r:     # r is a dict list
                if a['id'].startswith('ENS'):
                    ensg = (a['id'])
                    url = BASE_URL+ext_ebl+ensg+"?"
                    r = ebl_request(report, item, url, headers=headers)
                    if not r['seq_region_name'].startswith('CHR_'):
                        transcript = r['canonical_transcript'].split('.')[0]
                        symbol = r['display_name']
                        transcripts[transcript] = {'symbol':symbol, 'level': 'gene', 'given': item}
                        candidates_symbol.append(symbol)
            if len(candidates_symbol) > 1:
                report['multiple'].append({item: candidates_symbol})
    return transcripts


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
    for mesg in report['done']:
        print(f"  - {mesg}")

    if report['multiple']:
        print(f"{Color.BLUE}\n Multiple responses ({len(report['multiple'])}):")
        for mesg in report['multiple']:
            for k,v in mesg.items():
                print(f"  - {k}: {' '.join(v)}")

    if report['aborted']:
        print(f"{Color.PURPLE}\n Aborted ({len(report['aborted'])}):")
        for mesg in report['aborted']:
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
            fh.write(f"\n**Multiple Genes returned for one given ({len(report['multiple'])})**\n\n")
            for mesg in report['multiple']:
                for k,v in mesg.items():
                    fh.write(f"  - {k}: {' '.join(v)}")
        if report['aborted']:
            fh.write(f"\n\n**Genes/transcript missing ({len(report['aborted'])})**\n\n")
            for mesg in report['aborted']:
                fh.write(f"- {mesg}\n")


def gracefully_exit(args):
    pass


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
