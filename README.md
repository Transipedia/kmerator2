# Kmerator

## Prototype for decomposition of transcript or gene sequences and extraction of their specific k-mers

Kmerator is a prototype tool designed for the prediction of specific k-mers (also called tags) from input sequences, considering a reference genome and an ENSEMBL-like transcriptome. From these specific k-mers, it also outputs their corresponding specific contigs which are sequences of consecutive k-mers (overlapping length between k-mers must be k-1, otherwise, it's a new contig). Kmerator first uses Jellyfish [1] to create 2 requestable indexes from the reference genome and transcriptome, and second, decomposes your input transcript or gene sequences to count the occurences of each k-mer in the genome and transcriptome. Number of occurrences are then interpreted, in different manners, to select specific k-mer from your input. 

Kmerator strictly depends on a reference genome (fasta or jellyfish index format) and on an Ensembl fasta format transcriptome, you can find it there: https://www.ensembl.org/info/data/ftp/index.html. For a more complete k-mer filtering, we advice to merge the coding (cDNA) and non-coding (ncRNA) as one unique reference transcript.

Interpretation of occurrences counts differs depending on the ''level'' option : 

- gene 
- transcript 
- chimeras (must not be found in your reference genome and transcriptome)

If you want to use Kmerator on unannotated RNA signatures, you have to use the -u option. In this case, only "transcript" and "chimeras" levels are available.

## Dependencies

- Python >= v3.6
- Jellyfish >= 2.0


## Usage
```
kmerator.py [-h] (-s SELECTION [SELECTION ...] | -f FASTA_FILE) -g GENOME -t TRANSCRIPTOME   
			-l {gene,transcript,chimera} [-a APPRIS] [-u] [-k KMER_LENGTH] [--stringent]  
			[--threshold THRESHOLD] [-o OUTPUT] [-c CORES] [--verbose] [-v]
```

## arguments
```
optional arguments:
  -h, --help            show this help message and exit
  -s SELECTION [SELECTION ...], --selection SELECTION [SELECTION ...]
                        list of gene IDs or transcript IDs (ENSG, ENST or gene Symbol) to  
                        select inside your fasta transcriptome file and that you want to  
                        extract specific kmers from. If
                        you want to use your own sequences, you can give your fasta file   
                        with --fasta_file option.
  -f FASTA_FILE, --fasta-file FASTA_FILE
                        sequences in fasta format (genes or transcripts) that you want to  
                        extract specific kmers from. If you don't have your own sequences  
                        file, you can use a list of gene
                        IDs or transcript IDs with --selection option.
  -g GENOME, --genome GENOME
                        genome fasta file or jellyfish index (.jf) to use for k-mers   
                        request.
  -t TRANSCRIPTOME, --transcriptome TRANSCRIPTOME
                        transcriptome fasta file (ENSEMBL fasta format ONLY) to use for   
                        k-mers request and transcriptional variants informations.
  -l {gene,transcript,chimera}, --level {gene,transcript,chimera}
                        use 'gene', 'transcript' or 'chimera' to extract specific kmers at   
                        these different levels. Note that 'chimera' option activate the   
                        'unannotated' option.
  -a APPRIS, --appris APPRIS
                        indicate: 'homo_sapiens', 'mus_musculus', 'rattus_norvegicus',   
                        'danio_rerio', 'sus_scrofa', or virtually any specie available in   
                        APPRIS database [Rodriguez JM, et
                        al. Nucleic Acids Res. Database issue; 2017 Oct 23]. Genes have   
                        multiple possible transcripts called isoforms. This option selects   
                        the principal transcript defined
                        in APPRIS database. If this option is used and there is no data   
                        available, no connection or no principal transcript (non-coding   
                        genes), the longer transcript is
                        selected.
  -u, --unannotated     use this option if your provided sequences fasta file corresponds to     
  						annotations external from Ensembl. Otherwise, use ensembl   
  						annotations.
  -k KMER_LENGTH, --kmer-length KMER_LENGTH
                        k-mer length that you want to use (default 31).
  --stringent           FOR GENE LEVEL ONLY: use this option if you want to select
                        gene-specific k-mers present in ALL known transcripts for your gene.   
  					   If false, a k-mer is considered as gene-specific if present in at   
  					   least one isoform of your gene of interest.
  --threshold THRESHOLD
                        FOR GENE LEVEL ONLY: minimum fraction of annotated transcripts, for   
                        a given gene, containing this kmer to keep it (default: 0)
  -o OUTPUT, --output OUTPUT
                        output directory (default: 'output')
  -c CORES, --cores CORES
                        run n cores simultaneously (default: 1)
  --verbose             if you want some details while Kmerator is running.
  -v, --version         show program's version number and exit
```

## References

[1] Guillaume Marçais, Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers, Bioinformatics, Volume 27, Issue 6, 15 March 2011, Pages 764–770, https://doi.org/10.1093/bioinformatics/btr011
[2] Rodriguez JM, et al. Nucleic Acids Res. Database issue; 2017 Oct 23
