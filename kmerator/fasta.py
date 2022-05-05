### fasta.py


import sys


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
