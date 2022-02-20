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
    transcriptome_dict = load_transcriptome(args)
    ### build jellyfis genome and transcriptome
    jf_genome, jf_dir = run_jellyfish(args)

    ### get Ensembl name ID
    ensg_id = get_ensembl_id('RUNX1')

    ### get canonical transcript
    canonical_transcript = get_canonical_transcript(ensg_id)
    print(canonical_transcript)

    print(f"{Color.CYAN}\n     🪚  WORK IN PROGRESS. 🛠 next step : define 'find_longest_variant()' function.{Color.END}\n")

    print(f'''{Color.CYAN}
    Nota :
    - le dictionnaire renvoie t'il les mêmes positions que l'API d'ensembl ? (regarder avec gene_info.py)
    - Si oui, et en tenant compte du code complet, peut-être que le dictionnaire des transcripts ne sera pas utile...

    Sinon :
    - faire un git spécifique pour gene-info.py
    ''')



def load_transcriptome(args):
    """
    Load transcriptome file as dict
    if --unannotated is set -> basic conversion
    else description is modified
    """
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}"
                    f"Creating a dictionary of transcriptome file{Color.END}.")
    transcriptome = {}
    with open(args.transcriptome.name) as fh:
        seq = ""
        old_desc, new_desc = "", ""
        if args.unannotated:
            ### basic convertion to dict
            for line in fh:
                if line[0] == ">":
                    new_desc = line.rstrip()
                    if old_desc:
                        transcriptome[old_desc] = seq
                        seq = ""
                    old_desc = new_desc
                else:
                    seq += line.rstrip()
            transcriptome[old_desc] = seq
        else:
            ### convert fasta to dict with renamed description
            for line in fh:
                if line[0] == ">":
                    new_desc = line.split()
                    gene_name = new_desc[6].split(':')[1]
                    ensembl_transcript_name = new_desc[0].split('.')[0]
                    # ensembl_gene_name = new_desc[3].split(':')[1].split('.')[0]
                    new_desc = f"{gene_name}:{ensembl_transcript_name}"
                    if old_desc:
                        transcriptome[old_desc] = seq
                        seq = ""
                    old_desc = new_desc
                else:
                    seq += line.rstrip()
            transcriptome[old_desc] = seq
        if args.verbose: print(f"{Color.YELLOW}Transcript dictionary is done.{Color.END}")

    return transcriptome


def run_jellyfish(args):
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
        cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.cores}"
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



def get_ensembl_id(ref_gene):
    server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/"
    r = requests.get(server+ext+ref_gene+"?", headers={ "Content-Type" : "application/json"})
    if not r.ok:
        print(r.raise_for_status())
    if not r.json(): return None
    return [item['id'] for item in r.json() if item['id'].startswith('ENSG')][0]


def get_canonical_transcript(ensg_id):          ### APPRIS_function() in julia version
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/"
    r = requests.get(server+ext+ensg_id+"?", headers={ "Content-Type" : "application/json"})
    if not r.ok:
        return None
    return r.json()['canonical_transcript'].split('.')[0]


if __name__ == '__main__':
    main()
