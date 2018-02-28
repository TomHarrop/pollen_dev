#!/usr/bin/env python3

import os
import pathlib
import pandas
import shlex

#############
# FUNCTIONS #
#############

def resolve_path(my_path):
    return str(pathlib.Path(my_path).resolve())

###########
# GLOBALS #
###########

sample_key_file = 'data/sample_key.csv'
read_dir = 'data/fastq'
bbduk_adaptors = 'venv/bin/resources/adapters.fa'
bbduk_contaminants = 'venv/bin/resources/sequencing_artifacts.fa.gz'
star_reference_folder = 'output/010_ref/star_reference'


#########
# RULES #
#########

rule target:
    input:
        'output/010_ref/star_reference/Genome'


# 01 prepare STAR reference
rule star_reference:
    input:
        fasta = 'data/ref/TAIR10_Chr.all.fasta',
        gtf = ('output/010_ref/'
               'Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf')
    output:
        'output/010_ref/star_reference/Genome'
    params:
        genome_dir = star_reference_folder
    threads:
        30
    log:
        'output/logs/010_ref/star_reference.log'
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.genome_dir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 99 '
        '&> {log}'

rule preprocess_gtf:
    input:
        gtf = 'data/ref/Araport11_GFF3_genes_transposons.201606.gtf',
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff'
    output:
        gtf = ('output/010_ref/'
               'Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf')
    threads:
        1
    log:
        log = 'output/logs/010_ref/preprocess_gtf.log'
    script:
        'src/preprocess_gtf.R'


