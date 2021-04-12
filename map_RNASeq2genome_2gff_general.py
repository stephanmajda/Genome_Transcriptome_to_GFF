SAPLES_with_RNA_DATA =["199hm", "AR4D6","JBMS11", "LO226KS", "LO234KE","PR26KG","DS"]
from collections import defaultdict
import subprocess
from Bio import SeqIO
import re
import codecs

genome_folder = "/example/folder/"
trancriptome_folder = "../raw/transcriptome/"

outfolder ="../results/repeatscout/"

rule all:
    input:
        expand("{genome_folder}/{sample}.fa.masked", sample=SAPLES_with_RNA_DATA),
        expand("../results/augustus/{sample}_trinity_alignment.hint", sample=SAPLES_with_RNA_DATA),
        #expand("../results/augustus/{sample}_augustus_with_rnaseq_repeatmasked_cdna.fa", sample=SAPLES_with_RNA_DATA),
        expand("../results/minimap2/{sample}/{sample}_trinity_alignment.introns_added.gff", sample=SAPLES_with_RNA_DATA)

# include RNAseq:
rule find_transcripts:
    input:
        genome = expand("{genome_folder}{sample}.fa", sample=SAPLES_with_RNA_DATA)
    output:
        rna_files = "../results/trinity/file_list.txt"
    run:
        shell("find {trancriptome_folder} -name *trinity_transcripts.fasta > {output.rna_files}")

rule create_paf_files:
    input:
        rna_files = "../results/trinity/file_list.txt",
        genome = expand("{genome_folder}/{sample}.fa.masked", sample=SAPLES_with_RNA_DATA)
    output:
        paf = "../results/minimap2/{sample}/{sample}_trinity_alignment.paf"
    run:
        with open(input.rna_files, "r") as infile:
            for line in infile:
                if any(x in line for x in SAMPLES):
                    rna = str(line).strip()
                    shell("minimap2 {input.genome} {rna} -c -L -x splice -G 80K -t 16 > {output.paf}")

rule create_paf2bed:
    input:
        paf = "../results/minimap2/{sample}/{sample}_trinity_alignment.paf"
    output:
        bed = "../results/minimap2/{sample}/{sample}_trinity_alignment.bed"
    run:
        shell("paftools.js splice2bed {input.paf} > {output.bed}")


rule create_bed2gff:
    input:
        bed = "../results/minimap2/{sample}/{sample}_trinity_alignment.bed",
        #genome = "../results/ncbi_upload/Poteriospumella_lacustris_WGS_assembly_{sample}.trim.fa"
    output:
        gff2 = "../results/minimap2/{sample}/{sample}_trinity_alignment.gff",
    shell:
        "python bed_to_gff_converter.py {input.bed} {output.gff2}"

rule add_introns2gff:
    input:
        gff2 = "../results/minimap2/{sample}/{sample}_trinity_alignment.gff",
    output:
        gff_new = "../results/minimap2/{sample}/{sample}_trinity_alignment.introns_added.gff",
    run:
        with open(input.gff2,"r") as infile, open (output.gff_new,"w") as outfile:
            last_stop = 0
            for line in infile:
                if len(line.strip()) == 0 or line.startswith("#"):
                    continue
                data = line.split("\t")
                if data[2] == "mRNA":
                    mRNA_start = int(data[3])
                    mRNA_stop = int(data[4])
                    meta_info = data[5:]
                elif data[2] == "exon":
                    exon_start = int(data[3])
                    exon_stop = int(data[4])
                    if exon_start != mRNA_start or exon_stop != mRNA_stop:
                        if exon_start > mRNA_start:
                            intron = "\t".join(("\t".join(data[:2]),"intron",str(last_stop),str(exon_start-1),"\t".join(meta_info).replace("mRNA","intron")))
                            outfile.write(intron)
                        if exon_stop < mRNA_stop:
                            last_stop = exon_stop +1
                outfile.write(line)


rule gff2hint:
    input:
        gff2 = "../results/minimap2/{sample}/{sample}_trinity_alignment.introns_added.gff",
    output:
        hint = "../results/augustus/{sample}_trinity_alignment.hint",
    run:
        with open(input.gff2,"r") as infile, open (output.hint,"w") as outfile:
            for line in infile:
                if line.startswith("#") or len(line.strip()) == 0:
                    new_line = line
                else:
                    if line.split("\t")[2] == "mRNA":
                        continue
                    new_line = "\t".join(line.split("\t")[0:8]+["group=est_xyz; source=E\n"])

                outfile.write(new_line)

# include repeats:
rule repeatscout:
    input:
        genome = "{genome_folder}{sample}.fa"
    output:
        "../results/repeatscout/{sample}_repeats.txt"
    threads: 18
    run:
        cmd = {}
        output_lmer_frequency = outfolder+sample+"lmer.freq"
        cmd[0] = "build_lmer_table -sequence {} -freq {}".format(input.genome, output_lmer_frequency)
        cmd[1] = "RepeatScout -sequence {} -output {}_repeats.txt -freq {}".format(input.genome,outfolder+sample,output_lmer_frequency)
        for i in range(0,len(cmd)):
            print(i)
            subprocess.call(cmd[i], shell=True)

rule RepeatMasker:
    input:
        genome = "{genome_folder}{sample}.fa",
        repeats = "../results/repeatscout/{sample}_repeats.txt"
    output:
        "{genome_folder}/{sample}.fa.masked"
    threads: 18
    run:
        shell("RepeatMasker {input.genome} -lib {input.repeats} -pa {threads} -x")


rule create_paf_files_for_repeats:
    input:
        repeats_files = "../results/repeatscout/{sample}_repeats.txt",
        genome = "../results/new_spades/{sample}/scaffolds.fasta.filter"
    output:
        paf = "../results/repeatscout/{sample}_repeats_alignment.paf"
    run:
        shell("minimap2 {input.genome} {input.repeats_files} -c -L -G 80K -t 16 > {output.paf}")

rule create_paf2bed_for_repeats:
    input:
        paf = "../results/repeatscout/{sample}_repeats_alignment.paf"
    output:
        bed = "../results/repeatscout/{sample}_repeats_alignment.bed"
    run:
        shell("paftools.js splice2bed {input.paf} > {output.bed}")

rule mask_repeats:
    input:
        repeats_files = "../results/repeatscout/{sample}_repeats_alignment.bed",
        genome = "{genome_folder}/{sample}.fa"
    output:
        genome_masked = "{genome_folder}/{sample}.fa.masked"
    run:
        shell("bedtools maskfasta -soft -fi {input.genome} -bed {input.repeats_files} -fo {output.genome_masked}")

# final prediction:
rule gene_prediction:
    input:
        hint = "../results/augustus/{sample}_trinity_alignment.hint",
        genome_masked = "{genome_folder}/{sample}.fa.masked",
    output:
        aug_out = "../results/augustus/{sample}_augustus_with_rnaseq_repeatmasked.out",
    threads:3
    run:
        shell("augustus --species=arabidopsis {input.genome_masked} --hintsfile={input.hint} --gff3=on --softmasking=1 --singlestrand=true  --UTR=off --alternatives-from-evidence=false --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg > {output.aug_out}")
        #todo: remove comments from gff,,,,~/miniconda3/envs/augustus/config/extrinsic$ less extrinsic.M.RM.E.W.cfg
        # --AUGUSTUS_CONFIG_PATH=/home/sm/miniconda3/envs/augustus/config/
        #augustus_extrinsic.cfg

# extract nuc genes:
rule gene_extraction:
    input:
        augustus_file ="../results/augustus/{sample}_augustus_with_rnaseq_repeatmasked.out",
        genome = "{genome_folder}/{sample}.fa.masked",
    output:
        cdna_seq = "../results/augustus/{sample}_augustus_with_rnaseq_repeatmasked_cdna.fa",
    threads:3
    run:
        contig_dict={}
        gene_dict={}
        for seq_record in SeqIO.parse(input.genome, "fasta"):
            contig_dict[seq_record.id]=str(seq_record.seq)

        with open(input.augustus_file, "r",encoding='cp1252') as infile, open(output.cdna_seq, "w") as outfile:
            for line in infile:
                if not line.startswith("#"):
                    data = line.split("\t")
                    if data[2] == "transcript":
                        gene = data[-1].strip().split("Parent=")[1]
                        outfile.write("\n>{}_{}\n".format(wildcards.sample,gene))

                    elif data[2] == "CDS":
                        start = int(data[3]) - 1
                        end = int(data[4]) -1
                        outfile.write(contig_dict[data[0]][start:end])
