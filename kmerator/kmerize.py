### kmerize.py

import os
import subprocess
import multiprocessing

from fasta import fasta2dict


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
