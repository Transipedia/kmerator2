#!/usr/bin/env python3

"""
From gene, transcript or sequences, find specific kmers.
"""


import sys
import os
import subprocess

from utils import usage, checkup_args, Color


def main():
    """ Main function"""
    args = usage()
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Args:\n{args}{Color.END}")
    checkup_args(args)
    transcriptome_dict = load_transcriptome(args)
    # print(*list(transcriptome_dict.items())[20:22], sep='\n')
    jf_genome, jf_dir = run_jellyfish(args)
    ################################################################
    # TODO function run_jellyfish at line 310
    ################################################################


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
    ### create jellyfish PATH DIR
    jf_dir = f"{args.output}/jellyfish_indexes/{args.kmer_length}"
    os.makedirs(jf_dir, exist_ok=True)

    ### Compute jellyfish on transcriptome
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
    jf_transcriptome_dest = f"{jf_dir}/{'.'.join(os.path.basename(args.transcriptome.name).split('.')[:-1])}.jf"
    if os.path.exists(jf_transcriptome_dest):
        if args.verbose:
            print(f"{Color.YELLOW}{jf_transcriptome_dest} already exists, "
            f"keep it (manually remove to update it).{Color.END}")
    else:
        print(f"{Color.YELLOW}Compute Jellyfish on transcriptome{Color.END}")
        cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.cores}"
               f" -o {jf_transcriptome_dest} {args.transcriptome.name}")
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        except subprocess.CalledProcessError:
            sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                     f"{cmd}{Color.END}")

    ### Compute jellyfish on genome if genome is fasta file
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
    ext = args.genome.name.split('.')[-1]
    if ext == "fa" or ext == "fasta":
        jf_genome_dest = f"{jf_dir}/{'.'.join(os.path.basename(args.genome.name).split('.')[:-1])}.jf"
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

    ### Ending
    if args.verbose:
        print(f"{Color.YELLOW}Transcriptome kmer index output: {jf_transcriptome_dest}\n"
              f"Jellyfish done.{Color.END}")
    return None, None


if __name__ == '__main__':
    main()
