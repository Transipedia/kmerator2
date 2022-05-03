# usage.py

"""
Options and arguments to bo passed on to the program.
"""

import os
import sys
import argparse
import info


def usage():
    """
        Help function with argument parser.
    """
    parser = argparse.ArgumentParser(
        description=info.DOC,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    query_type = parser.add_mutually_exclusive_group(required=True)
    query_type.add_argument('-s', '--selection',
                            help=(
                                "list of gene IDs or transcript IDs (ENST, ENS or gene Symbol) to "
                                "select inside your fasta transcriptome file and that you want to "
                                "extract specific kmers from. For genes, kmerator search specific "
                                " kmers along the gene. For transcripts, it search specific kmers "
                                "to the transcript. You can also give a file with yours "
                                "genes/transcripts separated by space, tab or newline. "
                                "If you want to use your own unannotated sequences, you "
                                "must give your fasta file with --fasta_file option."
                                ),
                            nargs='+',
                            )
    query_type.add_argument('-f', '--fasta-file',
                            # ~ type=argparse.FileType('r'),
                            help=(
                                "Use this option when yours sequences are unannonated or provided "
                                "by a annotation file external from Ensembl. Otherwise, "
                                "use --selection option."
                                ),
                            )
    parser.add_argument('-g', '--genome',
                        # ~ type=argparse.FileType('r'),
                        help=(
                            "genome fasta file or jellyfish index (.jf) to use for k-mers requests."
                            ),
                        required=True,
                        )
    parser.add_argument('-t', '--transcriptome',
                        # ~ type=argparse.FileType('r'),
                        help=(
                            "transcriptome fasta file (ENSEMBL fasta format ONLY) to use for "
                            "k-mers request and transcriptional variants informations."
                            ),
                        required=True,
                        )
    parser.add_argument('-j', '--jellyfish-transcriptome',
                        # ~ type=argparse.FileType('r'),
                        help=(
                            "if your transcriptome (-t option) has already been converted by "
                            "jellyfish as a jf file, this avoids redoing the operation (be careful,"
                            "it must be the same transcriptome!)."
                            ),
                        )
    parser.add_argument('-c', '--specie',  ### replace -a --appris in julia version
                        help=(
                            "indicate a specie referenced in Ensembl, to help, follow the link "
                            "https://rest.ensembl.org/documentation/info/species. You can use "
                            "the 'name', the 'display_name' or any 'aliases'. For example human, "
                            "homo_sapiens or homsap are valid."
                            ),
                        default='human',
                        )
    parser.add_argument('-k', '--kmer-length',
                        type=int,
                        help="k-mer length that you want to use (default 31).",
                        default=31,
                        )
    parser.add_argument('--chimera',
                        action='store_true',
                        help="Only if with --fasta-file option.",
                        )
    parser.add_argument('--stringent',
                        action='store_true',
                        help=(
                            "FOR ANNOTATED GENE ONLY: use this option if you want to select "
                            "gene-specific k-mers present in ALL known transcripts for your gene. "
                            "If false, a k-mer is considered as gene-specific if present in at "
                            "least one isoform of your gene of interest."
                            ),
                        )
    parser.add_argument('--threshold',
                        type=float,
                        help=(
                            "FOR ANNOTATED GENE ONLY: minimum fraction of annotated transcripts, "
                            "for a given gene, containing this kmer to keep it (default: 0)"
                            ),
                        default=0,
                        )
    parser.add_argument('-o', '--output',
                        help="output directory (default: 'output')",
                        default='output',
                        )
    parser.add_argument('-p', '--procs',
                        type=int,
                        help="run n processes simultaneously (default: 1)",
                        default=1,
                        )
    parser.add_argument('--verbose',
                        action='store_true',
                        help="if you want some details while Kmerator is running.",
                        )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'{parser.prog} v{info.VERSION}',
                        )
    ### Go to 'usage()' without arguments or stdin
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    return args


def checkup_args(args):
    """ test genome and transcriptome files. """
    ### check if query fasta file is present
    if args.fasta_file and not os.path.isfile(args.fasta_file):
        sys.exit(f"{Color.RED}Error: {args.fasta_file!r} not found{Color.END}.")
    ### Check query fasta file extension
    if args.fasta_file and args.fasta_file.split('.')[-1] not in ("fa", "fasta"):
        sys.exit(f" {Color.RED}Error: {os.path.basename(args.fasta_file)} " \
                f"Does not appears to be a fasta file.{Color.END}")
    ### check if genome is present
    if not os.path.isfile(args.genome):
        sys.exit(f"{Color.RED}Error: {args.genome!r} not found{Color.END}.")
    ### Check genome fasta/jellyfish file extension
    if args.genome.split('.')[-1] not in ("fa", "fasta", "fna", "jf", "jl"):
        sys.exit(f" {Color.RED}Error: {os.path.basename(args.genome.name)}" \
                f"Does not appears to be a fasta or jelifish file.{Color.END}")
    ### check if transcriptome is present
    if args.transcriptome and not os.path.isfile(args.transcriptome):
        sys.exit(f"{Color.RED}Error: {args.transcriptome!r} not found{Color.END}.")
    ### Check genome fasta/jellyfish file extension
    ### Check transcriptome fasta file extension
    # ~ if args.transcriptome.name.split('.')[-1] not in ("fa", "fasta"):
        # ~ sys.exit(f" {Color.RED}Error: {os.path.basename(args.transcriptome.name)}" \
                # ~ f" does not appears to be a fasta file.{Color.END}")
    ### Check Gene Select option
    if args.selection:
        for gene in args.selection:
            if gene.startswith('ENS') and '.' in gene:
                sys.exit(f"{Color.RED}ENSEMBL annotations with version (point) like "
                        f"ENSTXXXXXXXX.1 is forbidden, just remove the '.1'.{Color.END}"
                        )
    ### Check if transcriptome is at ensembl format
    # ~ if not is_transcriptome_ensembl_file(args.transcriptome.readline()):
        # ~ sys.exit(f"{Color.RED}Error: {os.path.basename(args.transcriptome.name)} not in ENSEMBL fasta format, "
                # ~ f"use ENSEMBL transcriptome fasta file.{Color.END}")
    ## CANONICAL option works only with the gene annotated level
    # ~ if args.specie and (args.level != "gene" or args.unannotated):
        # ~ sys.exit(f"{Color.RED}Error: CANONICAL option works only with the gene annotated level{Color.END}")
    ## --chimera level works only with --fasta-file option
    if args.chimera and not args.fasta_file:
        sys.exit(f"{Color.RED}Error: '--chimera' needs '--fasta-file' option.{Color.END}")
    ### unanotated option want fasta_file option
    # ~ if args.unannotated and not args.fasta_file:
        # ~ sys.exit(f"{Color.RED}Error: with unannotated option, 'fasta-file' option is required.{Color.END}")
    ### if jellyfish of transcriptome is provided, check it
    if args.jellyfish_transcriptome and not args.jellyfish_transcriptome.split('.')[-1] == 'jf':
        sys.exit(f"{Color.RED}Error: {os.path.basename(args.jellyfish_transcriptome)!r} does not seem to be in jellyfish format.")


'''
def is_transcriptome_ensembl_file(first_row):
    """Checks if the line matches the Ensembl transcriptome file format"""
    line = first_row.split()
    try:
        if (not line[6].startswith('gene_symbol:') or
            not line[0].startswith('>ENST') or
            not line[3].startswith('gene:')):
            return False
    except:
        return False
    return True
'''

class Color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'
