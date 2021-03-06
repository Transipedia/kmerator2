#!/usr/bin/env python3

import sys
import os
import requests
import collections
from concurrent.futures import ThreadPoolExecutor
import time
from utils import Color
import json


def main():
    """ Function doc """
    ### some variables
    report = {'aborted': [], 'done': [], 'multiple': []}

    ### get transcript from Ensembl API
    print("GET ENSEMBL INFO")
    ensembl = Ensembl(args, report)
    transcripts = ensembl.get_ENST()

    print(f"Selection ({len(args.selection)})")
    print(*args.selection)

    print(f"\ntranscripts found ({len(transcripts)})")
    for k,v in transcripts.items():
        print(f"  {k}: {v}")


    print("\nREPORT")
    print(' done:', report['done'])
    print(' aborted:', report['aborted'])
    print(' multiple:', report['multiple'])


class Ensembl:
    """ Class doc """
    url = "https://rest.ensembl.org"
    nb_threads = 15                     # Because Ensembl limits requests to 15 per second'

    def __init__(self, args, transcript2gene_dict, report):
        """ Class initialiser """
        self.transcript2gene_dict = transcript2gene_dict
        self.report = report
        self.args = args
        self.headers = {"Content-Type" : "application/json"}   # header for the query


    def get_ENST(self):
        """ Function doc """
        transcripts = {}                # the dict to return
        urls = []                       # list of dicts where keys are 'item', 'level' and 'url'
        symbols_urls = []

        ### Define Ensembl URLs for each given item
        for item in self.args.selection:
            if item.upper().startswith('ENST'):
                urls.append({'level': 'transcript',
                             'item': item,
                             'url': f"{Ensembl.url}/lookup/id/{item}?",
                            })
            elif item.upper().startswith('ENSG'):
                urls.append({'level': 'gene',
                             'item': item,
                             'url': f"{Ensembl.url}/lookup/id/{item}?",
                            })
            else:
                symbols_urls.append({'item': item,
                                    'url': f"{Ensembl.url}/xrefs/symbol/{self.args.specie}/{item}?"})

        ### fetch ENSGxxx from each GENE SYMBOL
        with ThreadPoolExecutor(Ensembl.nb_threads) as pool:
            ensg_symbol_list = list(pool.map(self.ebl_request, symbols_urls))

        ### reduce symbol list to ENSG
        for symbol in ensg_symbol_list:
            symbol['response'] = [a['id'] for a in symbol['response'] if a['id'].startswith('ENSG')]

        ### Define symbol urls for urls list
        for symbol in ensg_symbol_list:
            if len(symbol['response']) == 0:
                # ~ self.report['aborted'].append(f"{symbol['item']!r} not found by Ensembl API.")
                continue
            for ensg in symbol['response']:
                item = symbol['item']
                urls.append({'level': 'gene',
                              'item': item,
                              'url': f"{Ensembl.url}/lookup/id/{ensg}?",
                            })

        ### fetch Ensembl ENSTxxx matching with given GENES, ENSGxxx and ENSTxxx
        find_xtime = collections.defaultdict(list)
        with ThreadPoolExecutor(Ensembl.nb_threads) as pool:
            enst_list = list(pool.map(self.ebl_request, urls))
        for dic in enst_list:
            response = dic['response']
            if not response:
                continue
            level = dic['level']
            item = dic['item']
            # ~ if  not 'display_name' in response:
                # ~ self.report['aborted'].append(f"{item!r} not found by Ensembl API.")
                # ~ continue
            if response['seq_region_name'].startswith('CHR_'):
                if not dic['item'] in find_xtime:  # to avoid repetition when multiples ENST are found for once gene
                    self.report['aborted'].append(f"{item!r} found in an alternative chromosome (Ensembl API).")
                continue
            if level == 'transcript':
                transcript = response['id'].upper()
            else:
                transcript = response['canonical_transcript'].split('.')[0]
            ### create transcript dict
            if transcript in self.transcript2gene_dict:
                symbol = self.transcript2gene_dict[transcript]
                transcripts[transcript] = {'given': item, 'symbol': symbol, 'level': level}
                if level == 'gene':
                    find_xtime[item].append(symbol)
            else:
                self.report['aborted'].append(f"{transcript!r} not found in transcriptome ({item}).")
                # ~ transcripts[transcript] = {'given': item, 'level': level}

        ### add multiple ENSTxxx found per gene
        for item, names in find_xtime.items():
            if len(names) > 1:
                self.report['multiple'].append({item: names})

        return transcripts


    def ebl_request(self, query):
        start_time = time.time()                # Because Ensembl limits requests to 15 per second
        ### part of query is returned in response
        response = query
        try:
            r = requests.get(query['url'], headers=self.headers)
        except requests.ConnectionError as err:
            sys.exit(f"{Color.RED}\n Error: Ensembl is not accessible or not responding.{Color.END}")
        ### remove useless info
        response.pop('url')
        if not r.ok:
            # ~ self.report['aborted'].append(f"{query['item']!r}: something went wrong whis the request to Ensembl.")
            r = []
        else:
            r = r.json()
        if not r:
            self.report['aborted'].append(f"{query['item']!r} not found by the Ensembl API.")
            r = []
        elif 'error' in r:
            self.report['aborted'].append(f"{query['item']!r} not found by Ensembl API.")
            r = []
        ### Because Ensembl limits requests to 15 per second
        limit = max(0, 1-(time.time()-start_time)+0.05)
        time.sleep(limit)
        response['response'] = r

        return response


class args:
    """ Class doc """
    specie = 'human'
    # ~ selection = [ 'GLOP', 'BRCA2', 'NPM1', 'CHI3L1', 'VPS29', 'ENST00000621131', 'ENSG00000159216',
        # ~ 'ENST00000614774', 'hstf1', 'braf', 'tp53', 'kras', 'egfr', 'ret', 'met', 'alk',
        # ~ 'BRCA1', 'stk11', 'ros1', 'notch1', 'CD46', 'GNG11', 'ADAMTS4', 'PCDH9', 'ZNF703', 'ZNF462'
        # ~ 'ZNF569', 'SNORA23', 'TXNRD1', 'ERRFI1', 'MARVELD2', 'DEPTOR', 'ELN', 'ADAMTS7', 'NOTCH3'
    # ~ ]
    selection = ['hstf1', 'glOp', 'enst00000614774', 'enst00000621131', 'enst00006211431', 'npm1', 'ensg00000159216']
    # ~ selection = ['hstf1']


if __name__ == "__main__":
    main()
