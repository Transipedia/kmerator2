#!/usr/bin/env python3


import requests, sys

server = "https://rest.ensembl.org"
ext = "/lookup/id/"
genes = (
    'ENSG00000157764',
    'ENSG00000099899',
    'ENSG00000181163',
    'ENSG00000122025',
    'ENSG00000159216',
)

for gene in genes:
    r = requests.get(server+ext+gene+"?", headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    print(f"---\n{r.json()}")
