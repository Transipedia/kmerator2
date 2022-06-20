### sequences.py


import os
import sys
from concurrent.futures import ThreadPoolExecutor

from utils import Color
from fasta import fasta2dict


class Sequences:
    """ Class doc """

    def __init__(self, args, report, transcripts, transcriptome_dict=None):
        """ Class initialiser """
        self.args = args
        self.report = report
        self.transcripts = transcripts
        self.transcriptome_dict = transcriptome_dict
        self.removed_transcripts = [] # list of transcript to remove when not found in transcriptome


    def build(self):
        """
        create files for each transcript
        Notice:
            transcripts == best_transcripts  when genes/transcripts are known
            transcripts == unannotated_transcripts  for unannotated sequences
        """
        output_seq_dir = os.path.join(self.args.output, 'sequences')
        ### Whith --selection option
        if self.args.selection:
            ### Abort if no transcripts found
            if not self.transcripts:
                sys.exit(f"{Color.RED}\n Error: no sequence found for {self.args.selection}{Color.END}")
            ### create output directory structure
            os.makedirs(output_seq_dir, exist_ok=True)

            ### create fasta files for each gene/transcript of selection found in transcriptome
            transcripts_list = [ {**{'enst':k}, **v} for k,v in self.transcripts.items()]
            with ThreadPoolExecutor(self.args.procs) as pool:
                list(pool.map(self.mk_annot_seq, transcripts_list))
            for trscrpt in self.removed_transcripts:
                self.transcripts.pop(trscrpt)
        ### With --fasta-file option
        else:
            ### read fasta file
            if self.args.debug: print(f"{Color.YELLOW}{'-'*12}\n\nBuild sequences without transcriptome.\n{Color.END}")
            fastafile_dict = fasta2dict(self.args.fasta_file)
            ### Abort if dict empy
            if not fastafile_dict:
                sys.exit(f"{Color.RED}Error: no sequence found for {args.fasta_file}")
            ### create output directory structure
            os.makedirs(output_seq_dir, exist_ok=True)
            ### create fasta files for unannotated sequences
            fastafile_list = [ {'desc': k, 'seq':v} for k,v in fastafile_dict.items()]
            with ThreadPoolExecutor(self.args.procs) as pool:
                list(pool.map(self.mk_unannot_seq, fastafile_list))


    def mk_annot_seq(self, transcript):
        """
        find canonical transcript in transcriptome and write his sequence in file.
        """
        # ~ if not 'symbol' in transcript:
            # ~ print(f"{Color.RED}No symbol name found ({transcript}.{Color.END}")
        desc = f"{transcript['symbol']}:{transcript['enst']}"

        if desc in self.transcriptome_dict:
            seq = self.transcriptome_dict[desc]
            if len(seq) < self.args.kmer_length:
                self.report['warning'].append(f"{desc} has a sequence length < {self.args.kmer_length} => ignored")
                return
            ### create fasta files
            outfile = f"{transcript['symbol'].replace('.','_')}.{transcript['enst']}.fa"[:255].replace(' ', '_').replace('/', '@SLASH@')
            outfile = f"{self.args.output}/sequences/{outfile}"
            with open(outfile, 'w') as fh:
                fh.write(f">{transcript['symbol']}:{transcript['enst']}\n{seq}")
            ### check if file is present
            if not os.path.isfile(outfile):
                print(f"{Color.RED}Error: unable to create {outfile}.")
        ### When transcript is not found
        else:
            self.report['aborted'].append(f"{transcript['enst']!r} not found in provided transcriptome (gene: {transcript['symbol']})")
            self.removed_transcripts.append(transcript['enst'])


    def mk_unannot_seq(self, entry):
        """ Function doc """
        outfile = f"{entry['desc'].replace(' ', '_').replace('/', '@SLASH@')}.fa"[:255]
        outfile = os.path.join(self.args.output, 'sequences', outfile)
        if len(entry['seq']) < self.args.kmer_length:
            self.report['aborted'].append(f"{entry['desc']} sequence length < {self.args.kmer_length} => ignored.")
            return
        self.transcripts.append(entry['desc'])
        with open(outfile, 'w') as fh:
            fh.write(f">{entry['desc'][:79]}\n{entry['seq']}")
        ### check if file is present
        if not os.path.isfile(outfile):
            print(f"{Color.RED}Error: unable to create {outfile}.")
