# Kmerator2

## Prototype for decomposition of transcript or gene sequences and extraction of their specific k-mers


Kmerator is a prototype tool designed for the prediction of specific k-mers (also called tags) from input sequences, considering a reference genome and an ENSEMBL-like transcriptome. From these specific k-mers, it also outputs their corresponding specific contigs which are sequences of consecutive k-mers (overlapping length between k-mers must be k-1, otherwise, it's a new contig). Kmerator first uses Jellyfish [1] to create 2 requestable indexes from the reference genome and transcriptome, and second, decomposes your input transcript or gene sequences to count the occurences of each k-mer in the genome and transcriptome. Number of occurrences are then interpreted, in different manners, to select specific k-mer from your input. 

Kmerator strictly depends on a reference genome (fasta or jellyfish index format) and on an Ensembl fasta format transcriptome, you can find it there: https://www.ensembl.org/info/data/ftp/index.html. For a more complete k-mer filtering, we advice to merge the coding (cDNA) and non-coding (ncRNA) as one unique reference transcript.

Kmerator2 is new version of kmerator, written in python3 (julia with first version), several options have changed. It is compatible with last versions of Ensembl transcriptome (version 103 max for kmerator). The output is improved and a file report is produced. 


#### Specific kmers

![](img/specific-kmers.png)

#### Specific config

![](img/specific-contigs.png)

## Dependencies

- Python >= v3.6
- Jellyfish >= 2.0


## Installation

### Solution 1 (preferred)

Install with pip

```
pip3 install kmerator2
```

If installed as user, ensure the directory `$HOME/.local/bin` is in your $PATH.


### Solution 2

Installation from github

```
git clone https://github.com/Transipedia/kmerator2.git
cp kmerator2/kmerator/kmerator.py /usr/local/bin/kmerator2  # or somewhere in your $PATH
cp kmerator2/kmerator/ktools.py /usr/local/bin/ktools       # or somewhere in your $PATH
```


## Usage
```
kmerator2 [-h] (-s SELECTION [SELECTION ...] | -f FASTA_FILE) -g GENOME -t TRANSCRIPTOME   
			-l {gene,transcript,chimera} [-a APPRIS] [-u] [-k KMER_LENGTH] [--stringent]  
			[--threshold THRESHOLD] [-o OUTPUT] [-c CORES] [--verbose] [-v]
```


## How use kmerator2

There are two main cases:

- you find for specific k-mers for annotated genes or transcripts : use the `--selection` option, followed by:
	- the list of gene and/or transcripts
	- or a file with the list of genes/transcripts
- you find for specific k-mers of unannotated sequences : use the `--fasta-file` option, followed by a fasta file containing yours requests. In case of you focuses on chimeras, add the `--chimera` option

### Differences between genes and transcripts

- When you find for a gene (symbol or Ensembl name), kmerator fetch sequence of its canonical transcript, extracts kmers and keep those that found only in the gene.
- When you find for a transcript, kmerator only keeps the kmer found in the transcript, and only in that transcript. If isoforms completely cover the transcript, no kmer will be kept.

## arguments

```
Choose one of the two options:
  -s SELECTION [SELECTION ...], --selection SELECTION [SELECTION ...]
                        list of gene IDs or transcript IDs (ENST, ENS or gene Symbol) to
                        select inside your fasta transcriptome file and that you want to
                        extract specific kmers from. For genes, kmerator search specific kmers
                        along the gene. For transcripts, it search specific kmers to the
                        transcript. You can also give a file with yours genes/transcripts
                        separated by space, tab or newline. If you want to use your own
                        unannotated sequences, you must give your fasta file with --fasta_file
                        option.
  -f FASTA_FILE, --fasta-file FASTA_FILE
                        Use this option when yours sequences are unannonated or provided by a
                        annotation file external from Ensembl. Otherwise, use --selection
                        option.
Mandatory:
  -g GENOME, --genome GENOME
                        genome fasta file or jellyfish index (.jf) to use for k-mers requests.
  -t TRANSCRIPTOME, --transcriptome TRANSCRIPTOME
                        transcriptome fasta file (ENSEMBL fasta format ONLY) to use for k-mers
                        request and transcriptional variants informations.
optional arguments:
  -j JELLYFISH_TRANSCRIPTOME, --jellyfish-transcriptome JELLYFISH_TRANSCRIPTOME
                        if your transcriptome (-t option) has already been converted by
                        jellyfish as a jf file, this avoids redoing the operation (be
                        careful,it must be the same transcriptome!).
  -S SPECIE, --specie SPECIE
                        indicate a specie referenced in Ensembl, to help, follow the link
                        https://rest.ensembl.org/documentation/info/species. You can use the
                        'name', the 'display_name' or any 'aliases'. For example human,
                        homo_sapiens or homsap are valid.
  -k KMER_LENGTH, --kmer-length KMER_LENGTH
                        k-mer length that you want to use (default 31).
  --chimera             Only if with --fasta-file option.
  --stringent           FOR ANNOTATED GENE ONLY: use this option if you want to select gene-
                        specific k-mers present in ALL known transcripts for your gene. If
                        false, a k-mer is considered as gene-specific if present in at least
                        one isoform of your gene of interest.
  -o OUTPUT, --output OUTPUT
                        output directory (default: 'output')
  -p PROCS, --procs PROCS
                        run n process simultaneously (default: 1)
  -d, --debug           if you want some details while Kmerator is running.
  --keep                keep intermediate files (sequences, indexes, separate tags and contigs
                        files).
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```


**Nota**: `kmerator2` has lost the `--level` option specifying the level (gene, transcript or chimera). It's now semi-automatic: if you give a `gene symbol` or a `ENSGxxx` the level is **gene**, if you give a `ENSTxxx` the level is **transcript**. You can mix gene symbols, ENSGxxxx and ENSTxxx. You need to use the `--fasta-file` in association with `--chimera` to operate at `chimera` level.

## Examples

**Get specific kmers from a set of known gene and/or transcripts:**
```
kmerator2 --selection BRCA2  npm1 ENSG00000159216 ENST00000614774 \
	--transcriptome /indexes/Homo_sapiens.GRCh38.cdna+ncrna-altchr.fa \
	-g /indexes/GRCh38_with_MT.jf \
	-o my-output \
	-p 8
```
Notes
- Genes and transcripts are mixed, but keep in mind that the behaviour is different between the two: in the case of genes, kmerator looks for gene-specific kmers. For transcripts, only kmers specific to the transcript are retained, excluding kmers found in isoforms of the transcript.
- The list of genes/transcript could be in file, separated by spaces, tab or newlines, like `--selection my-gene-file`.
- kmerator will search the file `/indexes/Homo_sapiens.GRCh38.cdna+ncrna-altchr.jf`, If not exists, it creates it in the output.

**Get specific kmers of unannotated sequences/transcripts:**
```
kmerator2 --fasta-file my-unannotated-seqs.fa \
	--transcriptome /indexes/Homo_sapiens.GRCh38.cdna+ncrna-altchr.fa \
	-g /indexes/GRCh38_with_MT.jf \
	-o my-output \
	-p 8
```



## ktools, a companion tool for kmerator

ktools can help you with some kmerator related tasks. For example, build the transcriptome could be tricky and repetitive (updated quaterly).

```
usage: ktools [-h] [-v] {mk-transcripts} ...

positional arguments:
  {mk-transcripts}
    mk-transcripts  make transcriptome

optional arguments:
  -h, --help        show this help message and exit
  -v, --version     show program's version number and exit
```

### build transcriptome

ktools get **cDNA** and **ncRNA** Ensembl transcriptome fasta files (last release by default), it concatene this files and remove alternative chromosomes.




## References

[1] Guillaume Marçais, Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers, Bioinformatics, Volume 27, Issue 6, 15 March 2011, Pages 764–770, https://doi.org/10.1093/bioinformatics/btr011
[2] Rodriguez JM, et al. Nucleic Acids Res. Database issue; 2017 Oct 23
