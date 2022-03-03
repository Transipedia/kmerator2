#!/usr/bin/env python3

"""
From genes, transcripts or sequences, find specific kmers.
"""


import sys
import os
import subprocess
import requests

import info
from utils import usage, checkup_args, Color


def main():
    """ Main function"""
    args = usage()
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Args:\n{args}{Color.END}")
    checkup_args(args)
    ### load transcriptome as dict
    if args.unannotated:
        transcriptome_dict = fasta_as_dict(args, args.fasta_file)
    else:
        transcriptome_dict = ensembl_fasta_as_dict(args, args.transcriptome.name)
    ### build jellyfish genome and transcriptome
    jf_genome, jf_dir = _run_jellyfish(args)

    ''' For testing (TO DELETE)
    ### get Ensembl name ID
    for gene_name in args.selection:
        ens_id = get_ensembl_id(args, gene_name)
        print(f"{ens_id = }")
        ### get canonical transcript
        canonical_transcript = get_canonical_transcript(ens_id) if ens_id else None
        print(f"{canonical_transcript = }")
        ### get longest transcript
        find_longest_variant(args, gene_name, transcriptome_dict)
    '''

    ### build sequences
    build_sequences(args, transcriptome_dict)

    print(f"{Color.CYAN}\n     🪚  WORK IN PROGRESS. 🛠 current step : define 'build sequences()'{Color.END}")

    print(f'''{Color.CYAN}
     Sinon : faire un git spécifique pour gene-info.py
    ''')


def ensembl_fasta_as_dict(args, fastafile):
    """
    Load Ensembl fasta file as dict,
    It changes header as SYMBOL:ENSTxxx:ENSGxxx (without version)
    """
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}"
                    f"Creating a dictionary from {fastafile} (Ensembl fasta file).{Color.END}.")
    ensembl_fasta_dict = {}
    with open(fastafile) as fh:
        seq = ""
        old_desc, new_desc = "", ""
        ### convert fasta to dict with renamed description
        for line in fh:
            if line[0] == ">":
                new_desc = line.split()
                gene_name = new_desc[6].split(':')[1]
                ensembl_transcript_name = new_desc[0].split('.')[0].lstrip('>')
                ensembl_gene_name = new_desc[3].split(':')[1].split('.')[0]
                # ensembl_gene_name = new_desc[3].split(':')[1].split('.')[0]
                new_desc = f"{gene_name}:{ensembl_transcript_name}:{ensembl_gene_name}"
                if old_desc:
                    ensembl_fasta_dict[old_desc] = seq
                    seq = ""
                old_desc = new_desc
            else:
                seq += line.rstrip()
        ensembl_fasta_dict[old_desc] = seq
        if args.verbose: print(f"{Color.YELLOW}Ensembl dictionary from {fastafile} is done.{Color.END}")
    return ensembl_fasta_dict


def fasta_as_dict(args, fastafile):
    '''
    Basic convertion of fasta file to dict
    It keeps all the header
    '''
    fasta_dict = {}
    with open(fastafile.name) as fh:
        seq = ""
        old_desc, new_desc = "", ""
        for line in fh:
            if line[0] == ">":
                new_desc = line.rstrip().lstrip('>')
                if old_desc:
                    fasta_dict[old_desc] = seq
                    seq = ""
                old_desc = new_desc
            else:
                seq += line.rstrip()
    fasta_dict[old_desc] = seq
    return fasta_dict


def _run_jellyfish(args):
    genome = args.genome.name
    ### create jellyfish PATH DIR
    jf_dir = f"{args.output}/jellyfish_indexes/{args.kmer_length}"
    os.makedirs(jf_dir, exist_ok=True)

    ### Compute jellyfish on TRANSCRIPTOME
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
    jf_transcriptome_dest = f"{jf_dir}/{'.'.join(os.path.basename(args.transcriptome.name).split('.')[:-1])}.jf"
    if os.path.exists(jf_transcriptome_dest):
        if args.verbose: print(f"{Color.YELLOW}{jf_transcriptome_dest} already exists, remove it manually to update.{Color.END}")
    else:
        print(f"{Color.YELLOW}Compute Jellyfish on transcriptome{Color.END}")
        cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.procs}"
               f" -o {jf_transcriptome_dest} {args.transcriptome.name}")
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        except subprocess.CalledProcessError:
            sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                     f"{cmd}{Color.END}")

    ### Compute jellyfish on GENOME if genome is fasta file
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the genome.{Color.END}")
    ext = args.genome.name.split('.')[-1]
    if ext == "fa" or ext == "fasta":
        jf_genome = '.'.join(os.path.basename(args.genome.name).split('.')[:-1]) + '.jf'
        jf_genome_dest = os.path.join(jf_dir, jf_genome)
        if os.path.exists(jf_genome_dest):
            if args.verbose:
                print(f"{Color.YELLOW}{jf_genome_dest} already exists, "
                f"keep it (manually remove to update it).{Color.END}")
        else:
            print(f"{Color.YELLOW}Compute Jellyfish on genome{Color.END}.")
            cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.cores}"
                f" -o {jf_genome_dest} {args.genome.name}")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
            except subprocess.CalledProcessError:
                sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                        f"{cmd}{Color.END}")
    else:
        if args.verbose: print(f"{Color.YELLOW}Jellyfish genome index already provided.{Color.END}")
        jf_genome = os.path.basename(genome)
        jf_dir = os.path.dirname(genome)

    ### Ending
    if args.verbose:
        print(f"{Color.YELLOW}Transcriptome kmer index output: {jf_transcriptome_dest}\n"
              f"Jellyfish done.{Color.END}")
    return jf_genome, jf_dir


def _get_ensembl_id(args, ref_gene):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/{args.canonical}/"
    r = requests.get(server+ext+ref_gene+"?", headers={ "Content-Type" : "application/json"})
    if not r.ok: return None
    if not r.json(): return None
    return [item['id'] for item in r.json() if item['id'].startswith('ENS')][0]


def _get_canonical_transcript(ens_id):          ### APPRIS_function() in julia version
    server = "https://rest.ensembL.org"
    ext = "/lookup/id/"
    try:
        r = requests.get(server+ext+ens_id+"?", headers={ "Content-Type" : "application/json"})
    except requests.exceptions.ConnectionError:
        print(f"Warning: {server!r} not responding")
        return None
    if not r.ok: return None
    if not r.json(): return None
    return r.json()['canonical_transcript'].split('.')[0]


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


def _store_fasta_seq(args, gene_name, ensembl_transcript_name, seq):
    '''Print fasta sequence'''
    gene_name = gene_name.replace('/', '@SLASH@')
    outfile = os.path.join(args.output, 'sequences', f'{gene_name}.{ensembl_transcript_name}.fa')
    with open(outfile, 'w') as fh:
        fh.write(f">{gene_name}:{ensembl_transcript_name}\n{seq}")


def build_sequences(args, transcriptome_dict):
    '''
    build fasta files with sequences of genes/transcripts
    when --fasta-file args is set --> sequences are fetch from the fasta-file arg
    when --selection args is set  --> sequences are fetch from the transcriptome arg
    case 1: unanotated is not set
        subcase 1: level is transcript (--selection must be set)
            get sequences -from transcriptome)
        subcase 2: level is gene
            sub2case 1: level is TRANSCRIPT
            sub2case 2: level is GENE
                --fasta-file is set OR gene and transcript in --selection
                sub3case 1: canonical is set
                    get canonical transcript
                sub3case 2: canonical is not set
    case 2: unanotated is set

    '''
    ### Creating individual sequence files from input fasta file or genes/transcripts list
    if args.fasta_file:
        fastafile = ensembl_fasta_as_dict(args, args.fasta_file.name)
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nSplitting sequences of fasta input file.\n{Color.END}")
    else:
        fastafile = transcriptome_dict
        if args.verbose:  print(f"{Color.YELLOW}{'-'*12}\n\nCreating sequences from transcripts/genes list.\n{Color.END}")
    ### create output directory structure
    output_seq_dir = os.path.join(args.output, 'sequences')
    os.makedirs(output_seq_dir, exist_ok=True)

    ### when 'unanotated' option is set
    if args.unannotated:
        ### without ensembl annotations
        for desc,seq in fastafile.items():
            outfile = f"{desc.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            if os.path.isfile(os.join(output_seq_dir, outfile)):
                print(f"{Color.PURPLE}Warning: file {outfile!r} already exist => ignored")
            if len(seq) >= args.kmer_length:
                with open(os.path.join(output_seq_dir, outfile), 'w') as fh:
                    fh.write(f">{desc[:79]}\n{seq}")
            else:
                print(f"{Color.PURPLE}Warning: {gene_name!r} sequence length < {args.kmer_length} => ignored")



    ### when 'unanotated' option is not set
    else:
        ### with ensembl annotations
        genes_already_processed = set()
        genes_analysed = set()

        ### Remember fastafile could be '--fasta_file' option or 'transcriptome'
        for desc,seq in fastafile.items():
            gene_name, ensembl_transcript_name, ensembl_gene_name = desc.split(':')

            ### Transcript level
            if args.level == "transcript":
                ## transcript is in selection list
                if args.selection and ensembl_transcript_name in args.selection:
                    if len(seq) < args.kmer_length:
                        print(f"{Color.CYAN}{ensembl_transcript_name}: sequence length < {arg.kmer_length} => ignored{Color.END}")
                        continue
                    gene_name = gene_name.replace('/', '@SLASH@')
                    outfile = os.path.join(output_seq_dir, f'{gene_name}.{ensembl_transcript_name}.fa')
                    with open(outfile, 'w') as fh:
                        fh.write(f">{gene_name}:{ensembl_transcript_name}\n{seq}")
                else:
                    continue

            ### gene level
            elif args.level == "gene":
                ### testing if gene has not been already processed
                ### if --fasta-file OR gene and transcript in --selection AND gene not in genes_already_processed
                if args.fasta_file or (gene_name in args.selection or ensembl_gene_name in args.selection) and gene_name not in genes_already_processed:
                ### --canonical option is set
                    if not gene_name in genes_analysed:
                        ### --canonical option is set
                        if args.canonical:
                            canonical_transcript = _get_canonical_transcript(ensembl_gene_name)
                            if canonical_transcript:
                                print(f"{'-'*12}\nUSE CANONICAL TANSCRIPT") # TO DELETE
                                print(f" GENE NAME: {gene_name}") # TO DELETE
                                print(f" DESC: {desc.split(':')[:2]}") # TO DELETE
                                print(f" CANONICAL TRANSCRIPT: {canonical_transcript}") # TO DELETE
                                print(f" ENSEMBL TRANSCRIPT NAME: {ensembl_transcript_name}") # TO DELETE
                                _store_fasta_seq(args, gene_name, canonical_transcript, seq)
                                genes_analysed.add(gene_name)
                                continue
                            else:
                                print(f"{Color.PURPLE}Warning: something went wrong with the Ensembl API for {gene_name!r}, using longest transcript.{Color.END}")
                        ### without canonical option, use longest transcript
                        longest_transcript = _find_longest_variant(args, gene_name, transcriptome_dict)
                        print(f"{'-'*12}\nUSE LONGEST TRANSCRIPT") # TO DELETE
                        print(f" GENE NAME: {gene_name}") # TO DELETE
                        print(f" LONGEST TRANSCRIPT: {longest_transcript}") # TO DELETE
                        if len(seq) >= args.kmer_length:
                            ### write results
                            _store_fasta_seq(args, gene_name, longest_transcript, seq)
                            genes_analysed.add(gene_name)
                        else:
                            print(f"{Color.PURPLE}Warning: {gene_name!r} sequence length < {args.kmer_length} => ignored")










if __name__ == '__main__':
    main()
