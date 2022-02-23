#!/usr/bin/env python3

import requests, sys

server = "https://rest.ensembl.org"
ext = "/xrefs/symbol/homo_sapiens/"
gene = 'BRCA2'


r = requests.get(server+ext+gene+"?", headers={ "Content-Type" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()

# ~ print(r.json())
for item in r.json():
    if item['type'] == 'gene' and item['id'].startswith('ENSG'):
        print(f"{gene}: {item['id']}")
