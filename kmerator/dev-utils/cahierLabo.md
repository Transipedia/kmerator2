# Décryptage de kmerator.jl


## Modifs

### Fonction checkup_variables()

Elle peut modifier des variables :

`fastafile` = fichier donné par --selection, sinon par --transcriptome 

J'ai déporté cette condition au début de build séquence.


### Fonction load_transcriptome()

* kmerator julia utilise la bibliothèqe FastIO.FastaReader pour charger les fasta, mais il est obligé de faire des modifs ensuite.
* kmerator python n'utilise pas de bibliothèque exterieure
* Le transcriptome peut être chargé 2 fois (lorsque l'otption `--fasta-file`est utilisée)

==> j'ai remplacé load_transcriptome() par les fonctions 
* `ensembl_fasta_as_dict()`charge le fichier fasta avec comme entête SYMBOL:ENSTxxx:ENSGxxx (sans n° de version)
* `fasta_as_dict()`. charge le fichier fasta en prenant la totalité des entêtes 



### Fonction APPRIS_function()

* lent
* Je l'ai changé par `get_canonical_transcript()` en utilisant l'API de Ensembl


## Déroulé



```
- args = usage()
- checkup_args()
- si option unanotated transcriptome_dict = fasta_as_dict()
  sinon transcriptome_dict = ensembl_fasta_as_dict()
- jf_genome, jf_dir = run_jellyfish
- build_séquence()
	- si fasta-file option -> fastafile = ensembl_fasta_as_dict()
	  sinon (selection option) --> fastafile = transcriptome_dict
	- si option unanotated
		pass
	  sinon (pas l'option unanotaded)
	  	- pour chaque ligne du dictionnaire fastafile on récupère desc (entete) et seq (séquence)
	  		- gene_name = gene symbol
	  		- ensembl_transcript_name = nom Ensembl du transcript
	  		- ensembl_gene_neme = nom Ensembl du gene
	  		
	  		- Si option level est à 'transcript'
	  			- si la longueur de la séquence < 31 (en fait args.length) --> on saute à la prochaine ligne
	  			- si 'ensembl_transcript_name' est dans l'option 'selection' --> on créé un fichier SYMBOL.ENSTxxx.fa contenant la séquence
	  		      sinon 'ensembl_transcript_name' n'est pas dans l'option 'selection' --> on saute à la prochaine ligne
	  			 
```

## Questions

