#!/usr/bin/env python3

"""
Give gene names, it returns general info on them.
"""

import sys
import os
import argparse
import requests


__appname__   = "gene-info"
__shortdesc__ = "Give gene names, it returns infos on them."
__licence__   = "none"
__version__   = "0.1.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    """ Function doc """
    args = usage()
    server = "https://rest.ensembl.org"
    for gene in args.genes:
        ### find type (ensembl, ncbi ..., otherwise symbol name)
        type = get_gene_type(gene)
        ### try to find info on gene
        info = fetch_info(args, gene, type)
        ### output results
        output(gene, type, info)


def output(gene, type, info):
    if info:
        print(f"Input: {gene} ({TYPES[type]})")
        ### for nice alignment
        lgth = 0
        for item in info:
            lgth = max(len(item), lgth)
        ### print for ensembl gene
        if type in ['ENSG', 'ENST', 'SYMBOL']:
            print(f"    {'Symbol Name':<{lgth}} {info['display_name']}")
            print(f"    {'Ensembl ID':<{lgth}} {info['id']}")
            print(f"    {'Type':<{lgth}} {info['object_type']}")
            print(f"    {'Specie':<{lgth}} {info['species']}")
            print(f"    {'Assembly':<{lgth}} {info['assembly_name']}")
            print(f"    {'Biotype':<{lgth}} {info['biotype']}")
            if 'description' in info: print(f"    {'Description':<{lgth}} {info['description']}")
            if 'canonical_transcript' in info: print(f"    {'Canonical Transcript':<{lgth}} {info['canonical_transcript']}")
            print(f"    {'Coordinates':<{lgth}} {info['seq_region_name']}:{info['start']}-{info['end']}")
            print(f"    {'Strand':<{lgth}} {info['strand']}")
            # ~ print(f"    {'Source':<{lgth}} {info['source']}")
            # ~ print(f"    {'Logic Name':<{lgth}} {info['logic_name']}")
            # ~ print(f"    {'Version':<{lgth}} {info['version']}")
            # ~ print(f"    {'DB type':<{lgth}} {info['db_type']}")
            # ~ if 'is_canonical' in info: print(f"    {'Canonical':<{lgth}} {info['is_canonical']}")
        else:
            for k,value in info.items():
                print(f"    {k:<{lgth}} {value}")
    else:
        print(f"Input: {gene} (NOT FOUND)")


def fetch_info(args, gene, type):
    '''
    Depending on the type (ensembl, ncbi, symbol name), try to find information about the gene
    '''
    ensg = gene
    if type == 'SYMBOL':
        ensg = _get_ENSG_from_SYMBOL(args, gene)
    info = _ensembl_lookup(ensg) if ensg else None
    return info


def _ensembl_lookup(gene):
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/"
    try:
        r = requests.get(server+ext+gene+"?", headers={ "Content-Type" : "application/json"})
    except:
        sys.exit('GLOP PAS')
    if not r.ok:
        return None
    return r.json()


def _get_ENSG_from_SYMBOL(args, gene):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/{args.specie}/"
    try:
        r = requests.get(server+ext+gene+"?", headers={ "Content-Type" : "application/json"})
    except:
        sys.exit('PAS GLOP')
    if 'error' in r.json():
        sys.exit(r.json()['error'])
    if not r.ok:
        sys.exit(r.raise_for_status())
    if not r.json(): return None
    ensg = [item['id'] for item in r.json() if item['id'].startswith('ENS')][0]
    return ensg


def get_gene_type(gene):
    '''
    return
    '''
    if gene.startswith('ENSG'): return 'ENSG'
    if gene.startswith('ENST'): return 'ENST'
    else: return 'SYMBOL'


TYPES = {
    'SYMBOL': 'Gene Symbol',        # Default
    ### Ensembl
    'ENSG': 'Ensembl Gene ID',
    'ENST': 'Ensembl Transcript ID',
    # ~ 'ENSP': 'Ensembl Protein ID',
    ### from http://idmap.genestimuli.org/
    # ~ 'N/A': 'NCBI gene ID',
    # ~ 'N/A': 'NCBI RefSeq ID',
    # ~ 'N/A': 'NCBI UniGene ID',
    # ~ 'N/A': 'Accession Number',
    # ~ 'N/A': 'UniProt ID',
    # ~ 'N/A': 'PDB ID',
    # ~ 'N/A': 'Prosite ID',
    # ~ 'N/A': 'PFam ID',
    # ~ 'N/A': 'InterPro ID',
    # ~ 'N/A': 'OMIM ID',
    # ~ 'N/A': 'PharmGKB ID',
    # ~ 'N/A': 'Affymetrix Probeset',
    # ~ 'N/A': 'HUGO Gene ID',
    ### RefSeq: NCBI Reference Sequence Database :
    ### https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
    # ~ 'AC_': '',
    # ~ 'NC_': '',
    # ~ 'NG_': '',
    # ~ 'NT_': '',
    # ~ 'NW_': '',
    # ~ 'NZ_': '',
    # ~ 'NM_': '',
    # ~ 'NR_': '',
    # ~ 'XM_': '',
    # ~ 'XR_': '',
    # ~ 'AP_': '',
    # ~ 'NP_': '',
    # ~ 'YP_': '',
    # ~ 'XP_': '',
    # ~ 'WP_': '',

}

def usage():
    """Help function with argument parser."""
    doc_sep = '=' * min(49, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    genes = '' if sys.stdin.isatty() else sys.stdin.read().rstrip().split()
    ### OPTION
    parser.add_argument("genes",
                        help="gene list or stdin (example: RUNX1 ENSG00000157764",
                        nargs='*' if genes else '+',
                        default=genes,
                        metavar=('gene'),
                       )
    parser.add_argument("-s", "--specie",
                        help="specie",
                        metavar="specie",
                        default='homo_sapiens',
                       )
    '''
    ### ARGUMENT WITHOUT OPTION
    parser.add_argument('--verbose',
                        action="store_true",          # boolean
                        help="Increase volubility",
                       )
    ### ARGUMENT WITH PREDIFINED OPTION
    parser.add_argument("-n", "--nombre",
                        type=int,
                        choices=[1, 2, 3],
                        help="un chiffre de 0 à 2",
                       )
    '''
    ### VERSIONNING
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1 and sys.stdin.isatty():
        sys.exit(parser.print_help())
    return parser.parse_args()

if __name__ == "__main__":
    main()
