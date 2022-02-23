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
    # ~ transcriptome_dict = load_transcriptome(args)
    if args.unannotated:
        transcriptome_dict = fasta_as_dict(args, args.transcriptome.name)
    else:
        transcriptome_dict = ensembl_fasta_as_dict(args, args.transcriptome.name)
    ### build jellyfis genome and transcriptome
    jf_genome, jf_dir = run_jellyfish(args)

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
     Sinon :
    - faire un git spécifique pour gene-info.py
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
    with open(fastafile) as fh:
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
                    new_desc = line.rstrip().lstrip('>')
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
                    ensembl_transcript_name = new_desc[0].split('.')[0].lstrip('>')
                    ensembl_gene_name = new_desc[3].split(':')[1].split('.')[0]
                    # ensembl_gene_name = new_desc[3].split(':')[1].split('.')[0]
                    new_desc = f"{gene_name}:{ensembl_transcript_name}:{ensembl_gene_name}"
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


def get_ensembl_id(args, ref_gene):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/{args.canonical}/"
    r = requests.get(server+ext+ref_gene+"?", headers={ "Content-Type" : "application/json"})
    if not r.ok: return None
    if not r.json(): return None
    return [item['id'] for item in r.json() if item['id'].startswith('ENS')][0]


def get_canonical_transcript(ens_id):          ### APPRIS_function() in julia version
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/"
    r = requests.get(server+ext+ens_id+"?", headers={ "Content-Type" : "application/json"})
    if not r.ok: return None
    if not r.json(): return None
    return r.json()['canonical_transcript'].split('.')[0]


def find_longest_variant(args, gene_name, transcriptome_dict):
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Finding the longest variant for the gene {gene_name}.{Color.END}")
    gene_name = gene_name.replace("/", "@SLASH@")           # sometimes, gene name contains a '/'
    # ~ for k,v in transcriptome_dict.items():
        # ~ print(k,v)
        # ~ sys.exit()
    variants_dict = { k:len(v) for (k,v) in transcriptome_dict.items() if k.startswith(f"{gene_name}:")}
    # ~ print(*[k for k in variants_dict], sep='\n')
    nb_variants = len(variants_dict)
    if args.verbose: print(f"{Color.YELLOW}Number of variants: {nb_variants}")
    longest_variant = None
    length = 0
    for k,v in variants_dict.items():
        if v > length:
            length = v
            longest_variant = ':'.join(k.split(':')[:2])
    print(f"{longest_variant = }")
    return longest_variant


def build_sequences(args, transcriptome_dict):
    ### Creating individual sequence files from input fasta file or genes/transcripts list
    if args.fasta_file:
        fastafile = ensembl_fasta_as_dict(args, args.fasta_file.name)
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nSplitting sequences of fasta input file.\n{Color.END}")
    else:
        fastafile = transcriptome_dict
        if args.verbose:  print(f"{Color.YELLOW}{'-'*12}\n\nCreating sequences from transcripts/genes list.\n{Color.END}")
    ### create output directory structure
    os.makedirs(os.path.join(args.output, 'sequences'), exist_ok=True)

    ### if 'unanotated' option is set
    if args.unannotated:
        ### without ensembl annotations
        pass
    ### if 'unanotated' option is not set
    else:
        ### with ensembl annotations
        '''
        As alternative, use ENSEMBL API to get ensembl-gene-name and ensemb-canonical-name ???
        '''
        genes_already_processed = []
        genes_analysed = []

    for desc,seq in fastafile.items():
            gene_name, ensembl_transcript_name, ensembl_gene_name = desc.split(':')

            ### Transcript level
            if args.level == "transcript":
                ## transcript not in selection list
                if args.selection and ensembl_transcript_name in args.selection:
                    if len(seq) < args.kmer_length:
                        print(f"{ensembl_transcript_name}: sequence length < {arg.kmer_length} => ignored")
                        continue
                    print("================")
                    print(f"{ensembl_transcript_name} est dans la selection")
                    print(f"GENE NAME: {gene_name}")
                    print(f"ENSEMB TRANSCRIPT NAME: {ensembl_transcript_name}")
                    print(f"ENSEMB GENE NAME: {ensembl_gene_name}")
                    print(f"SEQ LENGTH: {len(seq)}")
                    gene_name = gene_name.replace('/', '@SLASH@')
                    outfile = os.path.join(args.output, 'sequences', f'{gene_name}.{ensembl_transcript_name}.fa')
                    with open(outfile, 'w') as fh:
                        fh.write(f">{gene_name}:{ensembl_transcript_name}\n{seq}")

                else:
                    # ~ print(f"{ensembl_transcript_name} n'est pas dans la selection")
                    continue

            ### Transcript level
            else:
                pass
            '''
            ## Transcript level
            if level == "transcript"
                if !isempty(select_option) && !(ensembl_transcript_name in select_option)
                   continue
                else
                    if length("$seq") >= kmer_length
                        if verbose_option println("$ensembl_transcript_name: sequence length >= $kmer_length => continue") end
                        gene_name = replace(gene_name, "/" => "@SLASH@") # some gene names can contain slash characters that break the processus
                        FastaWriter("$output/sequences/$gene_name:$ensembl_transcript_name.fa") do fwsequence
                            write(fwsequence, [">$gene_name:$ensembl_transcript_name", "$seq"])
                        end
                    else println("$ensembl_transcript_name: sequence length < $kmer_length => ignored") end
                end
            '''


if __name__ == '__main__':
    main()
