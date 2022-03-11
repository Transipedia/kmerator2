#!/usr/bin/env python3

"""
From genes, transcripts or sequences, find specific kmers.
"""


import sys
import os
import subprocess
import requests
import multiprocessing
import time   # TO DELETE

import info
from utils import usage, checkup_args, Color


BASE_URL = "https://rest.ensembl.org"


def main():
    """ Main function"""
    ### Handle arguments
    args = usage()
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Args:\n{args}{Color.END}")
    checkup_args(args)

    transcriptome_dict = None
    canonical_transcripts = None
    ### when selection option is set
    if args.selection:
        ### get transcripts using Ensembl API
        print(f" ✨✨ Fetch some information from Ensembl API.")
        canonical_transcripts = _get_ensembl_transcripts(args)
        ### Load transcriptome as dict (necessary to build sequences and to found specific kmers
        print(f" ✨✨ Load transcriptome.")
        transcriptome_dict = ebl_fasta2dict(args.transcriptome)
        ### Build sequence using Ensembl API
        print(f" ✨✨ Build sequences.")
        build_sequences(args, canonical_transcripts, transcriptome_dict)
    else:
        print(f" ✨✨ Build sequences.")
        build_sequences(args)

    ### get specific kmer with multithreading
    print(f" ✨✨ Define specific kmers.")
    kmers = SpecificKmers(args, transcriptome_dict, canonical_transcripts)

    print(f"{Color.CYAN}\n     🪚  WORK IN PROGRESS. 🛠 next step : compute specific kmers'{Color.END}")
    print(f"{Color.CYAN}        Et aussi : faire un git spécifique pour gene-info.py")
    sys.exit()

    '''
    ### load transcriptome and fasta-file as dict
    transcriptome_dict = None
    fastafile_dict = None
    if args.selection:
        transcriptome_dict = ensembl_fasta2dict(args, args.transcriptome.name)
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nCreate transcriptome dictionnary.\n{Color.END}")
    if args.fasta_file:
        fastafile_dict = fasta2dict(args.fasta_file)
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nCreate input fasta file dictionnary.\n{Color.END}")

    ### build jellyfish genome and transcriptome
    jf_genome, jf_dir = _run_jellyfish(args)

    ### build sequences
    build_sequences(args, transcriptome_dict, fastafile_dict)
    '''


def build_sequences(args, transcripts=None, transcriptome_dict=None):
    ''''
    create files for each transcript
    '''
    output_seq_dir = os.path.join(args.output, 'sequences')
    ### Whith --selection option
    if args.selection:
        ### Abort if no transcripts found
        if not transcripts:
            sys.exit(f"{Color.RED}Error: no sequence found for {args.selection}")
        ### create output directory structure
        os.makedirs(output_seq_dir, exist_ok=True)
        ### Get the sequences and create files for each of them
        for transcript,gene in transcripts.items():
            if transcript in transcriptome_dict:
                seq = transcriptome_dict[transcript]
                if len(seq) < args.kmer_length:
                    print(f"{Color.YELLOW}Warning: {desc!r} sequence length < {args.kmer_length} => ignored{Color.END}")
                    continue
                ### create fasta files
                outfile = f"{gene[0].replace('.','_')}.{transcript}.fa"[:255].replace(' ', '_').replace('/', '@SLASH@')
                outfile = f"{args.output}/sequences/{outfile}"
                with open(outfile, 'w') as fh:
                    fh.write(f">{gene[0]}:{transcript}\n{seq}")
            ### if transcript not found
            else:
                print(f"{Color.YELLOW} Warning: {gene[0]}/{transcript} not found in provided transcriptome.{Color.END}")
            '''
            ### As alternative, fetch sequences with Ensembl API
            ext_ebl = f'/sequence/id/{transcript}?type=cdna;species={args.specie}'
            r = requests.get(BASE_URL+ext_ebl, headers={ "Content-Type": "text/plain"})
            seq = r.text
            '''
    ### Whith --fasta-file option
    else:
        ### read fasta file
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nBuild sequences without transcriptome.\n{Color.END}")
        fastafile_dict = fasta2dict(args.fasta_file)
        for desc,seq in fastafile_dict.items():
            outfile = f"{desc.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            if len(seq) < args.kmer_length:
                print(f"{Color.YELLOW}Warning: {desc!r} sequence length < {args.kmer_length} => ignored{Color.END}")
                continue
            with open(os.path.join(output_seq_dir, outfile), 'w') as fh:
                fh.write(f">{desc[:79]}\n{seq}")


def _get_ensembl_transcripts(args):
    ''''
    Works with --selection option,
    - get canonical transcript and symbol name if ENSG is provided
    - get symbol name if ENST is provided return dict as format {ENST: SYMBOL}
    - get canonical transcript if symbol name is provided
    return dict as format {ENST: SYMBOL}
    '''
    transcripts = {}                                # the dict to return
    ebl_motifs = ['ENSG', 'ENST']                   # Accepted Ensembl motifs
    ext_symbol = f"/xrefs/symbol/{args.specie}/"    # extension to get ENSG with SYMBOL
    ext_ebl = "/lookup/id/"                         # extension to get info with ENSG/ENST
    headers={ "Content-Type" : "application/json"}  # header for the query
    for item in args.selection:
        ### When ENSEMBL GENE NAME is provided, get canonical transcript
        if item[:4] == 'ENSG':
            url = BASE_URL+ext_ebl+item+"?"
            r = ebl_request(item, url, headers=headers)
            if not r: continue
            transcript = r['canonical_transcript'].split('.')[0]
            symbol = r['display_name']
            transcripts[transcript] = [symbol, 'gene']
        ### When ENST is provided, get symbol
        elif item[:4] == 'ENST':
            url = BASE_URL+ext_ebl+item+"?"
            r = ebl_request(item, url, headers=headers)
            if not r: continue
            if  not 'display_name' in r:
                print(f"{Color.RED}Error: display name of {item!r} not found, it will not be processed.{Color.END}")
                continue
            transcript = item.split('.')[0]
            symbol = r['display_name'].split('-')[0]
            # ~ parent = r.json()['Parent']
            # ~ r = requests.get(BASE_URL+ext_ebl+parent+"?", headers=headers)
            # ~ transcript = r.json()['canonical_transcript'].split('.')[0]
            # ~ symbol = r.json()['display_name']
            transcripts[transcript] = [symbol, 'transcript']
        ### In other cases, item is considered as NAME_SYMBOL
        else:
            url = BASE_URL+ext_symbol+item+"?"
            r = ebl_request(item, url, headers=headers)
            if not r: continue
            for a in r:
                if a['id'].startswith('ENSG'):
                    ensg = (a['id'])
                    url = BASE_URL+ext_ebl+ensg+"?"
                    r = ebl_request(item, url, headers=headers)
                    transcript = r['canonical_transcript'].split('.')[0]
                    symbol = r['display_name']
                    transcripts[transcript] = [symbol, 'gene']
    return transcripts


def ebl_request(item, url, headers):
    r = requests.get(url, headers=headers).json()
    if not r:
        print(f"{Color.RED}Error: {item!r} not found by Ensembl API, it will not be processed.{Color.END}")
        return None
    if 'error' in r:
        print(f"{Color.RED}Error: {r['error']}, it will not be processed.{Color.END}")
        return None
    return r


def fasta2dict(fasta_file):
    '''
    Basic convertion of fasta file to dict
    It keeps all the header
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
                new_desc = line.rstrip().lstrip('>')
                if old_desc:
                    fasta_dict[old_desc] = seq
                    seq = ""
                old_desc = new_desc
            else:
                seq += line.rstrip()
    fasta_dict[old_desc] = seq
    return fasta_dict


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
                new_desc = line.split('.')[0].lstrip('>')
                if old_desc:
                    fasta_dict[old_desc] = seq
                    seq = ""
                old_desc = new_desc
            else:
                seq += line.rstrip()
    fasta_dict[old_desc] = seq
    return fasta_dict


class SpecificKmers:
    """ Class doc """

    def __init__(self, args, transcriptome_dict, canonical_transcripts):
        """ Class initialiser """
        self.args = args
        self.transcriptome_dict = transcriptome_dict
        self.canonical_transcripts = canonical_transcripts
        ### compute Jellyfish on genome and transcriptome if not exists
        print('JELLYFISH')
        self.jellyfish()
        ### Sequences files to analyse
        seq_files = os.listdir(os.path.join(args.output, 'sequences'))
        with multiprocessing.Pool(processes=args.procs) as pool:
            results = pool.map(self.build_kmers, seq_files)
            # ~ results.wait()
            # ~ for res in results.get():
                # ~ print(res)

    def build_kmers(self, seq_file):
        """ Doc """
        ### Define output file names
        if self.args.selection:     # When '--selection' option is set
            gene_name, transcript_name = seq_file.split('.')[:2]
            tag_file = f"{gene_name}-{transcript_name}-specific_kmers.fa"
            contig_file = f"{gene_name}-{transcript_name}-specific_contigs.fa"
        else:    # When '--fasta-file' option is set
            gene_name = transcript_name = seq_file.split('.')[0]
            tag_file = f"{gene_name}-specific_kmers.fa"
            contig_file = f"{gene_name}-specific_contigs.fa"


        ### take the transcript sequence for jellyfish query
        sequence_fasta = fasta2dict(os.path.join(self.args.output,'sequences', seq_file))
        for desc, seq in sequence_fasta.items():
            sequence = seq
        print(desc)
        time.sleep(5)   # TO DELETE
        return seq_file


    def jellyfish(self):
        args = self.args
        ########################################################
        ### Jellyfish on genome and transcriptome
        ########################################################
        genome = args.genome
        ### To create jellyfish PATH DIR
        jf_dir = f"{args.output}/jellyfish_indexes/{args.kmer_length}"
        mk_jfdir = lambda x: os.makedirs(x, exist_ok=True)

        ### building kmercounts dictionary from jellyfish query on the genome

        ### Compute jellyfish on TRANSCRIPTOME
        if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
        if not args.jellyfish_transcriptome:
            args.jellyfish_transcriptome = f"{jf_dir}/{'.'.join(os.path.basename(args.transcriptome).split('.')[:-1])}.jf"
            print(" ✨✨ Compute Jellyfish on the transcriptome, please wait...")
            mk_jfdir(jf_dir)
            cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.procs}"
                   f" -o {args.jellyfish_transcriptome} {args.transcriptome}")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
            except subprocess.CalledProcessError:
                sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                         f"{cmd}{Color.END}")

        ### Compute jellyfish on GENOME if genome is fasta file
        if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the genome.{Color.END}")
        ext = args.genome.split('.')[-1]
        if ext == "fa" or ext == "fasta":
            print(" ✨✨ Compute Jellyfish on the genome, please wait...")
            mk_jfdir(jf_dir)
            jf_genome = '.'.join(os.path.basename(args.genome).split('.')[:-1]) + '.jf'
            args.jellyfish_genome = os.path.join(jf_dir, jf_genome)
            if os.path.exists(args.jellyfish_genome):
                if args.verbose:
                    print(f"{Color.YELLOW}{args.jellyfish_genome} already exists, "
                    f"keep it (manually remove to update it).{Color.END}")
            else:
                print(f"{Color.YELLOW}Compute Jellyfish on genome{Color.END}.")
                cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.procs}"
                    f" -o {args.jellyfish_genome} {args.genome}")
                try:
                    subprocess.run(cmd, shell=True, check=True, capture_output=True)
                except subprocess.CalledProcessError:
                    sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                            f"{cmd}{Color.END}")
        ### When jellyfish genome already exists
        else:
            if args.verbose: print(f"{Color.YELLOW}Jellyfish genome index already provided.{Color.END}")
            args.jellyfish_genome = genome

        ### Ending
        if args.verbose:
            print(f"{Color.YELLOW}Transcriptome kmer index output: {jf_transcriptome_dest}\n"
                  f"Jellyfish done.{Color.END}")


        '''
        ## building kmercounts dictionary from jellyfish query on the genome
        println("Jellyfish query -s $output/sequences/$splitted_fasta_files $jf_dir/$jf_genome")
        kmercounts_genome = read(`jellyfish query -s "$output/sequences/$splitted_fasta_files" "$jf_dir/$jf_genome"`, String)
        if verbose_option remotecall(replace_line, 1 , X, "jellyfish query $splitted_fasta_files on $jf_dir/$jf_genome finished") end
        kmercounts_genome = split(kmercounts_genome, "\n")[1:end-1]
        kmercounts_genome_dict = Dict()
        for mer in kmercounts_genome
            mer = split(mer)
            seq = mer[1]
            kmercounts_genome_dict["$seq"] = mer[2]
        end
        ## building kmercounts dictionary from jellyfish query on the transcriptome
        kmercounts_transcriptome = read(`jellyfish query -s "$output/sequences/$splitted_fasta_files" "$transcriptome"`, String)
        if verbose_option remotecall(replace_line, 1 , X, "Jellyfish query on $transcriptome done") end
        kmercounts_transcriptome_dict = Dict()
        kmercounts_transcriptome = split(kmercounts_transcriptome, "\n")[1:end-1]
        for mer in kmercounts_transcriptome
            mer = split(mer)
            seq = mer[1]
            kmercounts_transcriptome_dict["$seq"] = mer[2]
        end
        '''





""" _run_jellyfish
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
"""


""" _get_canonical_transcript
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
"""


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


""" _store_fasta_seq
def _store_fasta_seq(args, gene_name, ensembl_transcript_name, seq):
    '''Print fasta sequence'''
    gene_name = gene_name.replace('/', '@SLASH@')
    outfile = os.path.join(args.output, 'sequences', f'{gene_name}.{ensembl_transcript_name}.fa')
    with open(outfile, 'w') as fh:
        fh.write(f">{gene_name}:{ensembl_transcript_name}\n{seq}")
"""


""" ensembl_fasta2dict
def ensembl_fasta2dict(args, fastafile):
    '''
    Load Ensembl fasta file as dict,
    It changes header as SYMBOL:ENSTxxx:ENSGxxx (without version)
    '''
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
"""


""" build_sequences
def build_sequences(args, transcriptome_dict, fastafile_dict):
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
                sub3case 1: specie is set
                    get canonical transcript
                sub3case 2: specie is not set
    case 2: unanotated is set

    '''
    ### Creating individual sequence files from input fasta file or genes/transcripts list
    # ~ if args.fasta_file:
        # ~ fastafile = ensembl_fasta2dict(args, args.fasta_file.name)
        # ~ if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nSplitting sequences of fasta input file.\n{Color.END}")
    # ~ else:
        # ~ fastafile = transcriptome_dict
        # ~ if args.verbose:  print(f"{Color.YELLOW}{'-'*12}\n\nCreating sequences from transcripts/genes list.\n{Color.END}")
    ### create output directory structure
    output_seq_dir = os.path.join(args.output, 'sequences')
    os.makedirs(output_seq_dir, exist_ok=True)

    ### when '--fasta-file' option is set
    if args.fasta_file:
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nBuild sequences without transcriptome.\n{Color.END}")
        ### without ensembl annotations
        for desc,seq in fastafile_dict.items():
            outfile = f"{desc.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            if os.path.isfile(os.path.join(output_seq_dir, outfile)):
                print(f"{Color.PURPLE}Warning: file {outfile!r} already exist => ignored")
            if len(seq) >= args.kmer_length:
                with open(os.path.join(output_seq_dir, outfile), 'w') as fh:
                    fh.write(f">{desc[:79]}\n{seq}")
            else:
                print(f"{Color.PURPLE}Warning: {gene_name!r} sequence length < {args.kmer_length} => ignored")

    ### when selection option is set
    else:
        if args.verbose: print(f"{Color.YELLOW}{'-'*12}\n\nBuild sequences using transcriptome.\n{Color.END}")
        ### with ensembl annotations
        genes_already_processed = set()
        genes_analysed = set()

        ### Remember fastafile could be '--fasta_file' option or 'transcriptome'
        for desc,seq in transcriptome_dict.items():
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
                ### --specie option is set
                    if not gene_name in genes_analysed:
                        ### --specie option is set
                        if args.specie:
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
                        ### without specie option, use longest transcript
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
"""







if __name__ == '__main__':
    main()
