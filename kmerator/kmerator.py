#!/usr/bin/env python3

"""
From genes, transcripts or sequences, find specific kmers.

TODO
- [x] BUG du RUNX1 non trouvé dans le transcriptome v99 -> pas vraiment un bug.
- [ ] finaliser une première version pour pypi
- [ ] requêtes Ensembl threadées
- [ ] requêtes Ensembl avec un timeout
- [ ] requêtes Ensembl remplacées par le plus long transcript ???
- [ ] SpecificKmers dans un fichier séparé ???
- [ ] nettoyage des prints et autres bouts de codes commentés
- [>] Ajouter une option --report pour créer un rapport en markdown (date, auteur, commande complète, path, requêtes traitées, requêtes non traitées, etc.)
- [ ] pouvoir fournir un fichier à --selection
- [ ] pouvoir fournir un fichier de configuration
- [ ] le nom de l'option --fasta-file n'est pas terrible
- [✔] les genes/transcripts non traités devrait s'afficher à la fin
- [✔] {transcript: (gene, level)} n'est pas très clair  ->  {transcript: (gene, level, provided_by)} ou provided_by est 'Ensembl' ou 'longest'
- [ ] la création des séquences peut être multithreadée (une bonne idée ?)
- [ ] ajouter une option --keep pour garder les fichiers intermédiaires
- [ ] ajouter une option pour aggréger les fichiers de résultats
- [ ] ajouter dans le Pypi une commande pour générer automatiquement un transcriptome (en fonction de l'espèce se serait bien)
- [ ] ajouter dans le Pypi une commande qui renvoie la liste des espèces gérées par Ensembl
- [ ] faire un paquet DEB (avec pypi2deb ou py2deb)
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
    ### Handle arguments
    args = usage()
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Args:\n{args}{Color.END}")
    checkup_args(args)

    report = {'warning': [], 'done': []}
    transcriptome_dict = {}
    best_transcripts = {}                       # when genes/transcripts annotated (--selection)
    unannotated_transcripts = []                # when transcripts are unannotated (--fasta-file)
    ### when --selection option is set
    if args.selection:
        ### get transcripts using Ensembl API
        print(f" 🧬 Fetch some information from Ensembl API.")
        best_transcripts = _get_ensembl_transcripts(args, report)
        ### Load transcriptome as dict (necessary to build sequences and to found specific kmers
        print(f" 🧬 Load transcriptome.")
        transcriptome_dict = ebl_fasta2dict(args.transcriptome)
        ### Build sequence using provided transcriptome
        print(f" 🧬 Build sequences.")
        build_sequences(args, report, best_transcripts, transcriptome_dict)
    ### when --fasta-file option is set
    else:
        print(f" 🧬 Build sequences.")
        build_sequences(args, report, unannotated_transcripts)

    ### get specific kmer with multithreading
    print(f" 🧬 Extract specific kmers.")
    kmers = SpecificKmers(args, report, transcriptome_dict, best_transcripts, unannotated_transcripts)

    ### show warnings
    print(f"{Color.CYAN}\n Done ({len(report['done'])}):")
    for mesg in report['done']:
        print(f"  - {mesg}")

    print(f"{Color.PURPLE}\n Warning ({len(report['warning'])}):")
    for mesg in report['warning']:
        print(f"  - {mesg}")
    print(f"{Color.END}")


    print(f"{Color.CYAN}\n     🪚  Penser à faire un git spécifique pour gene-info.py.{Color.END}")
    sys.exit()


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
            desc = f"{values[0]}:{transcript}"
            if desc in transcriptome_dict:
                seq = transcriptome_dict[desc]
                if len(seq) < args.kmer_length:
                    # ~ print(f"{Color.YELLOW}Warning: {desc!r} sequence length < {args.kmer_length} => ignored{Color.END}")
                    report['warming'].append(f"{desc!r} sequence length < {args.kmer_length} => ignored")
                    continue
                ### create fasta files
                outfile = f"{values[0].replace('.','_')}.{transcript}.fa"[:255].replace(' ', '_').replace('/', '@SLASH@')
                outfile = f"{args.output}/sequences/{outfile}"
                with open(outfile, 'w') as fh:
                    fh.write(f">{values[0]}:{transcript}\n{seq}")
            ### When transcript is not found
            else:
                # ~ print(f"{Color.YELLOW}Warning: {values[0]}/{transcript} not found in provided transcriptome.{Color.END}")
                report['warning'].append(f"{values[0]}/{transcript} not found in provided transcriptome.")
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
        for desc,seq in fastafile_dict.items():
            outfile = f"{desc.replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
            if len(seq) < args.kmer_length:
                # ~ print(f"{Color.YELLOW}Warning: {desc!r} sequence length < {args.kmer_length} => ignored{Color.END}")
                report['warning'].append(f"{desc!r} sequence length < {args.kmer_length} => ignored")
                continue
            transcripts.append(desc)
            with open(os.path.join(output_seq_dir, outfile), 'w') as fh:
                fh.write(f">{desc[:79]}\n{seq}")


def _get_ensembl_transcripts(args, report):
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
            r = ebl_request(report, item, url, headers=headers)
            if not r: continue
            transcript = r['canonical_transcript'].split('.')[0]
            symbol = r['display_name']
            transcripts[transcript] = [symbol, 'gene']
        ### When ENST is provided, get symbol
        elif item[:4] == 'ENST':
            url = BASE_URL+ext_ebl+item+"?"
            r = ebl_request(report, item, url, headers=headers)
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
            r = ebl_request(report, item, url, headers=headers)
            if not r: continue
            for a in r:
                if a['id'].startswith('ENSG'):
                    ensg = (a['id'])
                    url = BASE_URL+ext_ebl+ensg+"?"
                    r = ebl_request(report, item, url, headers=headers)
                    transcript = r['canonical_transcript'].split('.')[0]
                    symbol = r['display_name']
                    transcripts[transcript] = [symbol, 'gene']
    return transcripts


def ebl_request(report, item, url, headers):
    r = requests.get(url, headers=headers)
    r = r.json()
    if not r:
        report['warning'].append(f"{item} not found by Ensembl API, it will not be processed.")
        return None
    if 'error' in r:
        report['warning'].append(f"{r[error]}, it will not be processed.")
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
                gene_name = line.split()[6].split(':')[1]
                transcript_name = line.split('.')[0].lstrip('>')
                new_desc = f"{gene_name}:{transcript_name}"
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
            transcript_name = transcript[0]         # ENST00000001
            gene_name = transcript[1][0]            # TP53
            level = transcript[1][1]                 # 'gene' or 'transcript'
            seq_file = f"{gene_name}.{transcript_name}.fa"
            ### Define all variants for a gene
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


        # ~ print(f"   🧬 compute {seq_file}")


        ### take the transcript sequence for jellyfish query
        sequence_fasta = fasta2dict(os.path.join(self.args.output,'sequences', seq_file))

        ### building kmercounts dictionary using jellyfish query on the genome
        cmd = (f"jellyfish query -s {os.path.join(self.seq_files_dir,seq_file)} {self.args.genome}")
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
            # ~ print(f"{Color.YELLOW}Warning: an error occured in jellyfish query command for {seq_file}:\n  {cmd}{Color.END}")
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

        # ~ print(fasta_kmer_list)      # debug
        # ~ print(fasta_contig_list)    # debug

        ## write tag files
        if fasta_kmer_list:
            tags_outdir = os.path.join(self.args.output, 'tags', str(self.args.kmer_length))
            os.makedirs(tags_outdir, exist_ok=True)
            with open(os.path.join(tags_outdir, tag_file), 'w') as fh:
                fh.write("\n".join(fasta_kmer_list) + '\n')
        if fasta_contig_list:
            contigs_outdir = os.path.join(self.args.output, 'contigs', str(self.args.kmer_length))
            os.makedirs(contigs_outdir, exist_ok=True)
            with open(os.path.join(contigs_outdir, contig_file), 'w') as fh:
                fh.write("\n".join(fasta_contig_list) + '\n')

        mesg = (f"{gene_name}:{transcript_name} (as {level}).")

        return mesg


    ### Jellyfish on genome and transcriptome
    def jellyfish(self):
        args = self.args
        genome = args.genome
        ### To create jellyfish PATH DIR
        jf_dir = f"{args.output}/jellyfish_indexes/{args.kmer_length}"
        mk_jfdir = lambda x: os.makedirs(x, exist_ok=True)

        ### building kmercounts dictionary from jellyfish query on the genome

        ### Compute jellyfish on TRANSCRIPTOME
        if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
        if not args.jellyfish_transcriptome:
            args.jellyfish_transcriptome = f"{jf_dir}/{'.'.join(os.path.basename(args.transcriptome).split('.')[:-1])}.jf"
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
        if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the genome.{Color.END}")
        ext = args.genome.split('.')[-1]
        if ext == "fa" or ext == "fasta":
            print(" 🧬 Compute Jellyfish on the genome, please wait...")
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
