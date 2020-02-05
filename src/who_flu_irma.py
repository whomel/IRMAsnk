#Example
#source ~/env/.mig/bin/activate
#snakemake --snakefile who_flu_irma.py -j 48 --config ifq=/Data-RAID5/iseq/raw/illumina/iSeqRun8/ out=out_iseq/iSeqRun8_IRMA/ -np

import subprocess, sys, os, glob
from os.path import join
#from pathlib import path
from os.path import basename
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

#############################################################################
#                                                                           #
# Description                                                               #
#                                                                           #
#############################################################################

#Pipeline designed to obtain consensus secuencies from raw fastq data generating denovo assemblies.
#The steps are: QualityTrim -> Denovo assembly -> Rename & Generate plots

# Input parameters  ------------------------------------------------------------------------
#

#
IFQ = config["ifq"]

workspace= config["out"]

#split correction step
threads=4

#NCBI databases
#blastDB="/Data-RAID5/iseq/databases/influenza_ncbi_2Dec19/influenza.fna"
#incomplete entries discarded DB
blastDB="/Data-RAID5/iseq/databases/influenza_ncbi_2Dec19_fixNames/influenza2.fna"
#Full DB
#blastDB="/Data-RAID5/iseq/databases/ncbi/nt"

segmentsFlu=["HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"]

trimmomatic="/home/migrau/apps/Trimmomatic-0.39"

## Functions -------------------------------------------------------------------

def fixNames(fafile):
    res=""
    segment=""
    listSeg=["HA","MP","NA","NP","NS","PA","PB1","PB2"]
    sampleName=fafile.split("/")[-2].split("_")[0]
    for index,record in enumerate(SeqIO.parse(fafile, 'fasta')):
        for g in listSeg:
            tmpDescr=record.description
            if g in tmpDescr:
                segment=g
        res+=">"+sampleName+"_"+segment+"\n"+str(record.seq)+"\n"
    return res

# Rules ------------------------------------------------------------------------
# 

#PAIRED READS
SAMPLES, PAIR= glob_wildcards(IFQ+"/{sample}_L001_{pair}_001.fastq.gz")

rule all:
    input:
        expand(workspace+'assemblies/{sample}/contigs.fasta', sample=SAMPLES),
        expand(workspace+"qualtrim/{sample}.R1.paired.fastq", sample=SAMPLES),
        expand(workspace+"qualtrim/{sample}.R2.paired.fastq", sample=SAMPLES),
        expand(workspace+'assemblies/rename/{sample}.fa', sample=SAMPLES),
        expand(workspace+'circos/{sample}-covX.txt', sample=SAMPLES),
        expand(workspace+"circos/{sample}-coverage.png", sample=SAMPLES),
        expand(workspace+"circos/png/{sample}-coverage.png", sample=SAMPLES),
        expand(workspace+"circos/tmp/{sample}_depth.coverage",sample=SAMPLES)

# rule preQuality:
#     input:
#         faR1=expand(IFQ+"{{sample}}.{pair}.fastq.gz", pair=["R1"]),
#         faR2=expand(IFQ+"{{sample}}.{pair}.fastq.gz", pair=["R2"])
#     output:
#         R1out=workspace+"qualtrim/preCheck/{sample}" 
#     shell:"""
#         fastqc -f fastq test1_S1_L001.R1.paired.fastq test1_S1_L001.R2.paired.fastq -o here/ 
#     """

#QUALITY FILTER
rule filter:
   input:
        faR1=expand(IFQ+"{{sample}}_L001_{pair}_001.fastq.gz", pair=["R1"]),
        faR2=expand(IFQ+"{{sample}}_L001_{pair}_001.fastq.gz", pair=["R2"])
   output:
        R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2out=workspace+"qualtrim/{sample}.R2.paired.fastq",
        R1out_unpaired=workspace+"qualtrim/{sample}.R1.unpaired.fastq",
        R2out_unpaired=workspace+"qualtrim/{sample}.R2.unpaired.fastq"
   params:
       trimmo=trimmomatic
   shell:"""
      
      java -jar {params.trimmo}/trimmomatic-0.39.jar PE -phred33 {input.faR1} {input.faR2} {output.R1out} {output.R1out_unpaired} {output.R2out} {output.R2out_unpaired} ILLUMINACLIP:{params.trimmo}/adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:100 HEADCROP:10 TRAILING:5
  """ 

#Denovo assembly using IRMA paired.
rule irma:
    input:
        R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2out=workspace+"qualtrim/{sample}.R2.paired.fastq"
    output:
        contigs=workspace+'assemblies/{sample}/contigs.fasta'
    params:
        folder=workspace+'assemblies/{sample}/',
        afolder=workspace+'assemblies/',
        tmpfolder='{sample}'
    shell:"""
        IRMA FLU {input.R1out} {input.R2out} {params.tmpfolder}
        mv $PWD/{params.tmpfolder} {params.afolder}
        if [ -s {params.folder}*HA*.fasta ]
        then
            cat {params.folder}*.fasta > {output.contigs}
        else
            touch {output.contigs}
        fi
    """

rule rename:
    input:
        fasta=workspace+'assemblies/{sample}/contigs.fasta'
    output:
        fasta=workspace+'assemblies/rename/{sample}.fa'
    run:
        newFasta = fixNames(input.fasta)
        with open(output.fasta,'w') as fo:
            fo.write(newFasta)

#Custom script (createGraphfiles_Full_17Nov.py) to create the circos depth plots.
rule subTypeAndDepthPlots:
    input:
        fasta=workspace+'assemblies/rename/{sample}.fa',
        R1=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2=workspace+"qualtrim/{sample}.R2.paired.fastq"
    output:
        covX=workspace+"circos/{sample}-covX.txt",
        png=workspace+"circos/{sample}-coverage.png",
        pngend=workspace+"circos/png/{sample}-coverage.png",
        tmpend=workspace+"circos/tmp/{sample}_depth.coverage"
    params:
        tmp=workspace+'assemblies/rename/{sample}-covX.txt',
        alltxt=workspace+'assemblies/rename/{sample}*.txt',
        circos=workspace+"circos/",
        circosTMP=workspace+"circos/tmp/",
        consensus=workspace+'assemblies/rename/{sample}',
        tmpdepth=workspace+"circos/{sample}_depth.coverage",
        tmppng=workspace+"circos/{sample}-coverage.png",
        fasta=workspace+'assemblies/{sample}/contigs.fasta',
        ws=workspace
    shell:"""
        #identify subtypes
        python ../tools/groupSubtypeIRMA.py {params.fasta} {params.ws}
        #depthPlots
        python /home/yi-mo/pipelines/miguel/paired/tools/createGraphfiles_Full_17Nov.py {input.fasta} {input.R1} {input.R2}
        mv {params.tmp} {output.covX}
        mv {params.alltxt} {params.circos}
        cp /home/migrau/pipelines/miguel/paired/tools/*conf {params.circos}
        if [ -s {params.tmpdepth} ] 
        then
            python /home/yi-mo/pipelines/miguel/paired/tools/generate_covplot.py {output.covX}
        else
            touch {params.tmpdepth}
            touch {params.tmppng}
            python /home/yi-mo/pipelines/miguel/paired/tools/generate_covplot.py {output.covX}
        fi
        mv {params.consensus}_depth.* {params.circosTMP}
        mv {params.consensus}_reftemp* {params.circosTMP}
        cp {output.png} {output.pngend}
    """
