#!/usr/bin/env python3

"""
From genes, transcripts or sequences, find specific kmers.
"""


import sys
import os
import subprocess

from utils import usage, checkup_args, Color


def main():
    """ Main function"""
    args = usage()
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Args:\n{args}{Color.END}")
    checkup_args(args)
    ### load transcriptome as dict
    transcriptome_dict = load_transcriptome(args)
    ### build jellyfis genome and transcriptome
    jf_genome, jf_dir = run_jellyfish(args)

    print(f"{Color.CYAN}\n  🪚  WORK IN PROGRESS. 🛠 curent step : define 'get_canonical_transcript()' function.{Color.END}\n")

    ### fetch canonical transcript
    # APPRIS_function in julia version
    get_canonical_transcript(args, 'ENSG00000099899')


def load_transcriptome(args):
    """
    Load transcriptome file as dict
    if --unannotated is set -> basic conversion
    else description is modified
    """
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}"
                    f"Creating a dictionary of transcriptome file{Color.END}.")
    transcriptome = {}
    with open(args.transcriptome.name) as fh:
        seq = ""
        old_desc, new_desc = "", ""
        if args.unannotated:
            ### basic convertion to dict
            for line in fh:
                if line[0] == ">":
                    new_desc = line.rstrip()
                    if old_desc:
                        transcriptome[old_desc] = seq
                        seq = ""
                    old_desc = new_desc
                else:
                    seq += line.rstrip()
            transcriptome[old_desc] = seq
        else:
            ### convert fasta to dict with renamed description
            for line in fh:
                if line[0] == ">":
                    new_desc = line.split()
                    gene_name = new_desc[6].split(':')[1]
                    ensembl_transcript_name = new_desc[0].split('.')[0]
                    # ensembl_gene_name = new_desc[3].split(':')[1].split('.')[0]
                    new_desc = f"{gene_name}:{ensembl_transcript_name}"
                    if old_desc:
                        transcriptome[old_desc] = seq
                        seq = ""
                    old_desc = new_desc
                else:
                    seq += line.rstrip()
            transcriptome[old_desc] = seq
        if args.verbose: print(f"{Color.YELLOW}Transcript dictionary is done.{Color.END}")

    return transcriptome


def run_jellyfish(args):
    genome = args.genome.name
    ### create jellyfish PATH DIR
    jf_dir = f"{args.output}/jellyfish_indexes/{args.kmer_length}"
    os.makedirs(jf_dir, exist_ok=True)

    ### Compute jellyfish on TRANSCRIPTOME
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the transcriptome.{Color.END}")
    jf_transcriptome_dest = f"{jf_dir}/{'.'.join(os.path.basename(args.transcriptome.name).split('.')[:-1])}.jf"
    if os.path.exists(jf_transcriptome_dest):
        if args.verbose: print(f"{Color.YELLOW}{jf_transcriptome_dest} already exists, remove it manually to update.{Color.END}")
    else:
        print(f"{Color.YELLOW}Compute Jellyfish on transcriptome{Color.END}")
        cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.cores}"
               f" -o {jf_transcriptome_dest} {args.transcriptome.name}")
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        except subprocess.CalledProcessError:
            sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                     f"{cmd}{Color.END}")

    ### Compute jellyfish on GENOME if genome is fasta file
    if args.verbose: print(f"{'-'*9}\n{Color.YELLOW}Compute Jellyfish on the genome.{Color.END}")
    ext = args.genome.name.split('.')[-1]
    if ext == "fa" or ext == "fasta":
        jf_genome = '.'.join(os.path.basename(args.genome.name).split('.')[:-1]) + '.jf'
        jf_genome_dest = os.path.join(jf_dir, jf_genome)
        if os.path.exists(jf_genome_dest):
            if args.verbose:
                print(f"{Color.YELLOW}{jf_genome_dest} already exists, "
                f"keep it (manually remove to update it).{Color.END}")
        else:
            print(f"{Color.YELLOW}Compute Jellyfish on genome{Color.END}.")
            cmd = (f"jellyfish count -m {args.kmer_length} -s 1000 -t {args.cores}"
                f" -o {jf_genome_dest} {args.genome.name}")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
            except subprocess.CalledProcessError:
                sys.exit(f"{Color.RED}An error occured in jellyfish command:\n"
                        f"{cmd}{Color.END}")
    else:
        if args.verbose: print(f"{Color.YELLOW}Jellyfish genome index already provided.{Color.END}")
        jf_genome = os.path.basename(genome)
        jf_dir = os.path.dirname(genome)

    ### Ending
    if args.verbose:
        print(f"{Color.YELLOW}Transcriptome kmer index output: {jf_transcriptome_dest}\n"
              f"Jellyfish done.{Color.END}")
    return jf_genome, jf_dir


### APPRIS_function() in julia version
def get_canonical_transcript(args, ref_gene):
    url = f"http://apprisws.bioinfo.cnio.es/rest/exporter/id/{args.appris}/{ref_gene}?methods=appris&format=json&sc=ensembl"
    print(url)
    pass

''' --> JULIA VERSION
function APPRIS_function(gene_ref)
    ! verbose_option ? print("\r") : println("\r ------------")
    print("Finding the principal isoform based on APPRIS database for $gene_ref... \n")
    ## Request to appris
    url = "http://apprisws.bioinfo.cnio.es/rest/exporter/id/$APPRIS_option/" *
            "$gene_ref?methods=appris&format=json&sc=ensembl"
    if verbose_option println("\nAPPRIS url: $url") end
    ## make_API_call(url)

    function http_req(url)
        try
            r = HTTP.request("GET", url; verbose=0)
            res = String(r.body)
            res = JSON.Parser.parse(res)
            return res
        catch err
            if verbose_option println("ERROR:\n $err") end
            res = "NODATA"
            return res
        end
    end

    transcripts = http_req(url)
    if verbose_option && isa(transcripts, Array)
        for transcript in transcripts
            println("First transcript found:")
            for (key,value) in transcript
                println("  $key: $value")
            end
            break
        end
    end
    ## Finding the best isoform, find Principals
    principals = []
    for (i,value) in enumerate(transcripts)
        try
            if occursin("PRINCIPAL:", transcripts[i]["reliability"])
                push!(principals, transcripts[i]["reliability"])
            end
        catch
        end
    end
    ## if no principals, return NODATA
    if isempty(principals)
        if verbose_option
            println("No principal isoform detected, this function will return " *
                        "'NODATA' and the longest transcript will be selected")
        end
        return(transcripts)
    end
    if verbose_option println("principals: $principals \n") end
    ## Select best Principal (minimum level)
    levels = []
    for i in 1:length(principals)
        push!(levels, split(principals[i],":")[2])
    end
    levels = map(x-> parse(Int, x), levels)
    level = minimum(levels)
    ## if multiple transcripts with better Principal, select the biggest
    selected_transcripts = hcat(map(x -> try x["transcript_id"]
                                         catch end, transcripts),
                                map(x -> try parse(Int, x["length_na"])
                                         catch end, transcripts),
                                map(x -> try x["reliability"]
                                         catch end, transcripts))
    selected_transcripts = selected_transcripts[map(x -> x != nothing, selected_transcripts[:, 3]), :]
    selected_transcripts = selected_transcripts[map(x -> occursin("PRINCIPAL:$level",x), selected_transcripts[:, 3]), :]
    max_length = maximum(selected_transcripts[:,2])
    selected_transcripts = selected_transcripts[map(x -> x == max_length, selected_transcripts[:, 2]), :]
    best_transcript = String( unique(selected_transcripts[:,1])[1])
    if verbose_option println("APPRIS result : $best_transcript \nAPPRIS function finished \n\r ------------ \n\n") end
    return(best_transcript)
end # end of APPRIS function
'''


if __name__ == '__main__':
    main()
