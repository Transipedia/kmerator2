#!/usr/bin/env python3

"""
From genes, transcripts or sequences, find specific kmers.

TODO
- [ ] ADD: Handle alternative to jellyfish, like kmc and kmtricks
- [ ] voir si des noms HUGO ou GENE SYMBOL commencent par ENS ==> sinon problème car il cherchera un ENS ENSEMBL au lieu de GENE SYMBOL (ligne 213 : elif item.startswith('ENS'):)
- [ ] BUG: Si je lui donne un ENSG, je ne retrouve pas le nom ensuite (car transformé en gene-SYMBOL) --> mettre son nom dans 'transcripts' (ligne 222 : transcripts[transcript] = [symbol, 'gene'])
- [ ] finaliser une première version pour pypi
- [ ] requêtes Ensembl threadées
- [ ] requêtes Ensembl avec un timeout
- [ ] EVAL: requêtes Ensembl remplacées par le plus long transcript (sous-entendu : y a t'il des genes qui n'ont pas de transcripts canoniques)???
- [ ] SpecificKmers dans un fichier séparé ???
- [ ] nettoyage des prints et autres bouts de codes commentés
- [ ] pouvoir fournir un fichier à --selection
- [ ] pouvoir fournir un fichier de configuration
- [ ] le nom de l'option --fasta-file n'est pas terrible
- [ ] la création des séquences peut être multithreadée (une bonne idée ?)
- [ ] ajouter une option --keep pour garder les fichiers intermédiaires (et sinon les supprimer) - ou le contraire --del-details
- [ ] ajouter une option pour aggréger les fichiers de résultats
- [ ] ajouter dans le Pypi une commande qui renvoie la liste des espèces gérées par Ensembl
- [ ] faire un paquet DEB (avec pypi2deb ou py2deb)
- [✔] ajouter dans le Pypi une commande pour générer automatiquement un transcriptome (en fonction de l'espèce ce serait bien)
- [✔] Ajouter une option --report pour créer un rapport en markdown (date, auteur, commande complète, path, requêtes traitées, requêtes non traitées, etc.)
- [✔] faire un test pour la souris
- [✔] kmerator.py : si le jellyfish du transcriptome est au même endroit (et avec le meme nom) que le .fa, le gérer automatiquement.
- [✔] Vérifier que la casse est gérées
- [✔] aujourdh'ui, plantage sur Ensembl (ddos sur les dns ?) avec un "ConnectionError" ==> message warning + se rabattre sur les plus longs transcripts
- [x] BUG du RUNX1 non trouvé dans le transcriptome v99 -> pas vraiment un bug.
- [✔] supprimer les sous-répertoire '31'
- [x] BUG: le transcript ENST00000621131 (VPS29) n'a pas été fait ! Même si aucun kmer trouvé, il faudrait créer un fichier vide... Mis quoi mettre dans le merged ? plus besoin car dans le rapport
- [✔] les genes/transcripts non traités devraient s'afficher à la fin
- [✔] {transcript: (gene, level)} n'est pas très clair  ->  {transcript: (gene, level, provided_by)} ou provided_by est 'Ensembl' ou 'longest'
"""

"""
New kmerator:
speed:
    - testing with 5 genes, new version: 13 sec, old version: 63 sec
ergonomic:
    - no need to specify --level, you can mix symbol name, ENST name and ENSG.
    - '--selection' option accept both list of genes or file (with the list of genes)
    - if not specified, the jellyfish trancriptome is automatically searched in:
        - the same directory of the fasta transcriptome (with same name but '.jf' extension)
        - the output directory of the jellyfish transcriptome (useful when re-played kmerator)
lisibility:
    - kmerator now show missing genes
installation:
    - no need julia, it's python, so, to install `pip3 install kmerator2`
output:
    - a report formated in markdown stores some information like user, date, command, etc.
    - merged tags and contigs are available (directly in the outpout directory)
"""


import sys
import os
import subprocess
import requests
import multiprocessing
import shutil
import getpass
from datetime import datetime
# ~ import time   # TO DELETE

import info
from utils import usage, checkup_args, Color


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
        best_transcripts = get_ensembl_transcripts(args, report)
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


class SpecificKmers:
    """ Class doc """

    def __init__(self, args, report, transcriptome_dict, best_transcripts, unannotated_transcripts):
        """ Class initialiser """
        self.args = args
        self.rev = rev = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                          'a':'t', 'c':'g', 'g':'c', 't':'a'}       # reverse base
        ### create a shared dict among multiple processes with Manager()
        ### (show https://ourpython.com/python/multiprocessing-how-do-i-share-a-dict-among-multiple-processes)
        manager = multiprocessing.Manager()
        self.transcriptome_dict = manager.dict(transcriptome_dict)
        self.best_transcripts = manager.dict(best_transcripts)
        ### compute Jellyfish on genome and transcriptome if not exists
        self.jellyfish()
        ### Sequences files to analyse
        self.seq_files_dir = os.path.join(self.args.output, 'sequences')
        ### launch workers

        transcripts = best_transcripts.items() if args.selection else unannotated_transcripts
        with multiprocessing.Pool(processes=self.args.procs) as pool:
            mesg = pool.map(self.worker, transcripts)
            report['done'] += mesg


    def worker(self, transcript):
        '''
        transcript is a dict when '--selection' is set, else it is a list
        '''
        ### Define some variables: gene_name, transcript_name, variants_dic and output file names
        fasta_kmer_list = []                # specific kmers list
        fasta_contig_list = []              # specific contigs list
        ## When '--selection' option is set
        if self.args.selection:
            transcript_name = transcript[0]          # ENST00000001
            gene_name = transcript[1]['symbol']      # tp53
            level = transcript[1]['level']           # 'gene' or 'transcript'
            given_name = transcript[1]['given']      # TP53 (gen/transcript given by user)
            seq_file = f"{gene_name}.{transcript_name}.fa"
            ### Define all variants for a gene
            # ~ for k,a in self.transcriptome_dict.items():
                # ~ print(k,a)
                # ~ break
            variants_dict = { k:v for k,v in self.transcriptome_dict.items() if k.startswith(gene_name) }
            nb_variants = len(variants_dict)
            # ~ print(f"{gene_name}: {variants_dict.keys()}")
            tag_file = f"{gene_name}-{transcript_name}-{level}-specific_kmers.fa"
            contig_file = f"{gene_name}-{transcript_name}-{level}-specific_contigs.fa"
        ## When '--chimera' option is set
        elif self.args.chimera:
            seq_file = f"{transcript.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            gene_name = transcript_name = transcript
            level = 'chimera'
            tag_file = f"{gene_name}-chimera-specific_kmers.fa"
            contig_file = f"{gene_name}-chimera-specific_contigs.fa"
        ## When '--fasta-file' option is set
        else:
            seq_file = f"{transcript.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            gene_name = transcript_name = transcript
            level = 'transcript'
            tag_file = f"{gene_name}-transcript-specific_kmers.fa"
            contig_file = f"{gene_name}-transcript-specific_contigs.fa"

        ### take the transcript sequence for jellyfish query
        sequence_fasta = fasta2dict(os.path.join(self.args.output,'sequences', seq_file))

        ### building kmercounts dictionary using jellyfish query on the genome
        cmd = (f"jellyfish query -s {os.path.join(self.seq_files_dir,seq_file)} {self.args.jellyfish_genome}")
        try:
            kmercounts_genome = subprocess.run(cmd, shell=True, check=True, capture_output=True).stdout.decode().rstrip().split('\n')
        except subprocess.CalledProcessError:
            sys.exit(f"{Color.RED}Error: an error occured in jellyfish query command for {seq_file}:\n  {cmd}{Color.END}")
            return None
        kmercounts_genome_dict = {}
        for mer in kmercounts_genome:
            seq, count = mer.split()
            kmercounts_genome_dict[seq] = int(count)

        ### building kmercounts dictionary using jellyfish query on the transcriptome
        cmd = (f"jellyfish query -s {os.path.join(self.seq_files_dir,seq_file)} {self.args.jellyfish_transcriptome}")
        try:
            kmercounts_transcriptome = subprocess.run(cmd, shell=True, check=True, capture_output=True).stdout.decode().rstrip().split('\n')
        except subprocess.CalledProcessError:
            self.report['warning'].append(f"an error occured in jellyfish query command for {seq_file}:\n  {cmd}")
            return None
        kmercounts_transcriptome_dict = {}
        for mer in kmercounts_transcriptome:
            seq, count = mer.split()
            kmercounts_transcriptome_dict[seq] = int(count)

        ### initialization of count variables
        i = 0       # kmer number
        j = 1       # contig number
        total_kmers = len(kmercounts_transcriptome_dict)
        if self.args.verbose: print(f"[{seq_file}]: Total kmers in kmercounts_transcriptome_dict= {total_kmers}")

        ## creating a new dictionary with kmers and their first position in our query sequence
        kmer_starts = {}
        kmer_placed = 0

        for mer in kmercounts_transcriptome_dict:
            ### get the first position of the kmer in the sequence
            kmer_placed += 1
            kmer_starts[mer] = next(iter(sequence_fasta.values())).index(mer)

        if self.args.verbose: print(f"[{seq_file}]: Total kmers found in sequence_fasta = {len(kmer_starts)}")
        ### rearrange kmer_starts as list of sorted tuple like (position, kmer)
        kmer_starts_sorted = sorted(list(zip(kmer_starts.values(), kmer_starts.keys())))  # array sorted by kmer position
        # ~ position_kmer_prev = first(kmer_starts_sorted[1])
        position_kmer_prev = kmer_starts_sorted[0][0]
        contig_seq = "" # initialize contig sequence
        ### for each kmer, get the count in both genome and transcriptome
        kmers_analysed = 0
        for tuple in kmer_starts_sorted:
            ### from the kmer/position sorted list, we extract sequence if specific (occurence ==1)
            mer = tuple[1]              # kmer sequence
            position_kmer = tuple[0]    # kmer position
            # ~ startt = time.time()
            kmers_analysed += 1
            per = round(kmers_analysed/total_kmers*100)     # to show percentage done ?

            if mer in kmercounts_genome_dict.keys():
                genome_count = kmercounts_genome_dict[mer]
            else:
                revcomp_mer = [self.rev[base] for base in mer]
                genome_count = kmercounts_genome_dict[revcomp_mer]
            transcriptome_count = kmercounts_transcriptome_dict[mer]

            ### Case of annotated genes/transcripts
            if level == 'gene':
                ## if the kmer is present/unique or does not exist (splicing?) on the genome
                if genome_count <= 1:
                    variants_containing_this_kmer = [k for k,v in variants_dict.items() if mer in v]
                    if self.args.stringent and transcriptome_count == nb_variants == len(variants_containing_this_kmer):
                        # kmers case
                        i += 1
                        tmp = len(variants_containing_this_kmer)
                        fasta_kmer_list.append(f">{gene_name}-{transcript_name}.kmer{i} ({tmp}/{nb_variants})\n{mer}")
                        # contigs case
                        if i == 1:
                            contig_seq = mer
                            position_kmer_prev = position_kmer
                        elif i>1 and position_kmer == position_kmer_prev+1:
                            contig_seq = f"{contig_seq}{mer[-1]}"
                            position_kmer_prev = position_kmer
                        else:
                            fasta_contig_list.append(f">{gene_name}-{transcript_name}.contig{j}\n{contig_seq}")
                            j = j+1
                            contig_seq = mer
                            position_kmer_prev = position_kmer
                    elif not self.args.stringent and transcriptome_count == len(variants_containing_this_kmer) and transcriptome_count > nb_variants * self.args.threshold:
                        ### kmers case
                        i += 1
                        tmp = len(variants_containing_this_kmer)
                        fasta_kmer_list.append(f">{gene_name}-{transcript_name}.kmer{i} ({tmp}/{nb_variants})\n{mer}")
                        ### contigs case
                        if i == 1:
                            contig_seq = mer
                            position_kmer_prev = position_kmer
                        elif i > 1 and position_kmer == position_kmer_prev + 1:
                            contig_seq = f"{contig_seq}{mer[-1]}"
                            position_kmer_prev = position_kmer
                        else:
                            fasta_contig_list.append(f">{gene_name}-{transcript_name}.contig{j}\n{contig_seq}")
                            j += 1
                            contig_seq = mer
                            position_kmer_prev = position_kmer

            elif level == 'transcript':
                if self.args.fasta_file and transcriptome_count == 0 and genome_count <= 1:
                    ### kmers case
                    i += 1
                    fasta_kmer_list.append(f">{gene_name}.kmer{i}\n{mer}")
                    ### contigs case
                    if i == 1:
                        contig_seq = mer
                        position_kmer_prev = position_kmer
                    elif i > 1 and position_kmer == position_kmer_prev+1:
                        contig_seq = f"{contig_seq}{mer[-1]}"
                        position_kmer_prev = position_kmer
                    else:
                        fasta_contig_list.append(f">{gene_name}.contig{j}\n{contig_seq}")
                        j += 1
                        contig_seq = mer
                        position_kmer_prev = position_kmer

                elif self.args.selection and transcriptome_count == 1 and genome_count <= 1:
                    ### kmers case
                    i += 1
                    fasta_kmer_list.append(f">{gene_name}-{transcript_name}.kmer{i}\n{mer}")
                    # ~ print(f">{gene_name}-{transcript_name}.kmer{i}\n{mer}")                     # TO DELETE
                    ### contigs case
                    if i == 1:
                        contig_seq = mer
                        position_kmer_prev = position_kmer
                    elif i > 1 and position_kmer == position_kmer_prev + 1:
                        contig_seq = f"{contig_seq}{mer[-1]}"
                        position_kmer_prev = position_kmer
                    else:
                        fasta_contig_list.append(f">{gene_name}.contig{j}\n{contig_seq}")
                        # ~ print(f">{gene_name}.contig{j}\n{contig_seq}")                              # TO DELETE
                        j += 1
                        contig_seq = mer
                        position_kmer_prev = position_kmer

            ### Case of unannotated sequences
            elif level == 'chimera':
                if transcriptome_count == genome_count == 0:
                    ### kmers case
                    i += 1
                    fasta_kmer_list.append(f">{gene_name}.kmer{i}\n{mer}")
                    ### contig case
                    if i == 1:
                        contig_seq = mer
                        position_kmer_prev = position_kmer
                    elif i > 1 and position_kmer == position_kmer_prev + 1:
                        contig_seq = f"{contig_seq}{mer[-1]}"
                        position_kmer_prev = position_kmer
                    else:
                        fasta_contig_list.append(f">{gene_name}.contig{j}\n{contig_seq}")
                        j += 1
                        contig_seq = mer
                        position_kmer_prev = position_kmer
            else:
                sys.exit(f"{Color.RED}Error: level {level} unknown")

        ### append last contig in list
        if level == "gene" and contig_seq:
            fasta_contig_list.append(f">{gene_name}-{transcript_name}.contig{j}\n{contig_seq}")
        elif level == "transcript" and self.args.selection and contig_seq:
            fasta_contig_list.append(f">{gene_name}-{transcript_name}.contig{j}\n{contig_seq}")
        elif (level == "chimera" or (level == "transcript" and self.args.fasta_file)) and contig_seq:
            fasta_contig_list.append(f">{gene_name}.contig{j}\n{contig_seq}")

        ## write tag files
        if fasta_kmer_list:
            tags_outdir = os.path.join(self.args.output, 'tags')
            os.makedirs(tags_outdir, exist_ok=True)
            with open(os.path.join(tags_outdir, tag_file), 'w') as fh:
                fh.write("\n".join(fasta_kmer_list) + '\n')
        ## write contig files
        if fasta_contig_list:
            contigs_outdir = os.path.join(self.args.output, 'contigs')
            os.makedirs(contigs_outdir, exist_ok=True)
            with open(os.path.join(contigs_outdir, contig_file), 'w') as fh:
                fh.write("\n".join(fasta_contig_list) + '\n')

        if self.args.selection:
            mesg = (f"{gene_name}:{transcript_name} as {level} level ({given_name}).")
        else:
            mesg = (f"{gene_name} as {level} level.")
        return mesg


    ### Jellyfish on genome and transcriptome
    def jellyfish(self):
        args = self.args
        genome = args.genome
        ### To create jellyfish PATH DIR
        jf_dir = f"{args.output}/jellyfish_indexes"
        mk_jfdir = lambda x: os.makedirs(x, exist_ok=True)

        ### building kmercounts dictionary from jellyfish query on the genome

        ### Compute jellyfish on TRANSCRIPTOME
        if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
        root_path = '.'.join(args.transcriptome.split('.')[:-1])
        root_basename = os.path.basename(root_path)
        jelly_candidate = f"{root_path}.jf"
        jelly_dest = f"{jf_dir}/{root_basename}.jf"
        ### check for existing jellyfish transcriptome
        if not args.jellyfish_transcriptome:
            ### at the same location of fasta transcriptome
            if os.path.isfile(jelly_candidate):
                args.jellyfish_transcriptome = jelly_candidate
            ### where the jellyfich file must be created
            if os.path.isfile(jelly_dest):
                args.jellyfish_transcriptome = jelly_dest
        ### do jellyfish on transcriptome fasta file
        if not args.jellyfish_transcriptome:
            tr_root_file = '.'.join(os.path.basename(args.transcriptome).split('.')[:-1])
            args.jellyfish_transcriptome = f"{jf_dir}/{tr_root_file}.jf"
            print(" 🧬 Compute Jellyfish on the transcriptome, please wait...")
            mk_jfdir(jf_dir)
            cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.procs}"
                   f" -o {args.jellyfish_transcriptome} {args.transcriptome}")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
            except subprocess.CalledProcessError:
                sys.exit(f"{Color.RED}An error occured in jellyfish count command:\n"
                         f"{cmd}{Color.END}")

        ### Compute jellyfish on GENOME if genome is fasta file
        ext = args.genome.split('.')[-1]
        if ext == "fa" or ext == "fasta":
            mk_jfdir(jf_dir)
            jf_genome = '.'.join(os.path.basename(args.genome).split('.')[:-1]) + '.jf'
            args.jellyfish_genome = os.path.join(jf_dir, jf_genome)
            if os.path.exists(args.jellyfish_genome):
                if args.verbose:
                    print(f"{Color.YELLOW}{args.jellyfish_genome} already exists, "
                    f"keep it (manually remove to update it).{Color.END}")
            else:
                print(" 🧬 Compute Jellyfish on the genome, please wait...")
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

    print(f"{Color.CYAN}\n     🪚  Penser à mettre gene-info.py dans bio2m-dev-laboratory.{Color.END}")


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





if __name__ == '__main__':
    main()
