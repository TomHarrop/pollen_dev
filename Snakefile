#!/usr/bin/env python3

import os
import pathlib
import pandas
import shlex


#############
# FUNCTIONS #
#############

def find_input_files(wildcards):
    my_stage = wildcards.stage
    my_plant = wildcards.plant
    my_samples = list(
        sample_key[(sample_key.stage == my_stage) &
                   (sample_key.plant == my_plant)]['filename'])
    my_r1 = [os.path.join(read_dir, x)
             for x in my_samples if '_R1_' in x][0]
    return {'r1': my_r1}


def resolve_path(my_path):
    return str(pathlib.Path(my_path).resolve())


###########
# GLOBALS #
###########

sample_key_file = 'data/sample_key_all.csv'
read_dir = 'data/fastq_all'
bbduk_adaptors = 'data/bbmap_resources/adapters.fa'
bbduk_contaminants = 'data/bbmap_resources/sequencing_artifacts.fa.gz'
star_reference_folder = 'output/010_ref/star_reference'

# containers
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
star_container = 'shub://TomHarrop/singularity-containers:star_2.6.0c'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
bio_container = 'shub://TomHarrop/singularity-containers:biopython_1.72'

#########
# SETUP #
#########

# generate name to filename dictionary
sample_key = pandas.read_csv(sample_key_file)

#########
# RULES #
#########

rule target:
    input:
        'output/090_deseq/dds.Rds'

# 09 DESeq analysis
rule generate_deseq_object:
    input:
        gene_calls = 'output/080_filter-background/gene_calls.csv',
        detected_genes = 'output/080_filter-background/detected_genes.csv'
    output:
        dds = 'output/090_deseq/dds.Rds'
    threads:
        1
    log:
        log = 'output/logs/090_deseq/generate_deseq_object.log'
    singularity:
        r_container
    script:
        'src/generate_deseq_object.R'

# 08 filter genes by expression
rule filter_backgroud:
    input:
        tpm = 'output/070_tpm/tpm.csv',
        intergenic_tpm = 'output/060_cutoffs/intergenic_tpm.csv'
    output:
        detected_genes = 'output/080_filter-background/detected_genes.csv',
        gene_calls = 'output/080_filter-background/gene_calls.csv',
        freqpoly = 'output/080_filter-background/cutoff_poly.pdf',
        violin = 'output/080_filter-background/cutoff_violin.pdf'
    log:
        log = 'output/logs/080_filter-background/filter_genes.log'
    threads:
        1
    singularity:
        r_container
    script:
        'src/filter_genes.R'


# 07 calculate TPM
rule calculate_tpm:
    input:
        count_files = expand(
            ('output/030_star-pass2/{stage}_{plant}.ReadsPerGene.out.tab'),
            stage=['UNM', 'PUNM', 'BCP', 'TCP'],
            plant=['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']),
        gtf = ('output/010_ref/'
               'Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf')
    params:
        star_dir = 'output/030_star-pass2'
    output:
        tpm = 'output/070_tpm/tpm.csv'
    threads:
        1
    log:
        log = 'output/logs/070_tpm/tpm.log'
    singularity:
        r_container
    script:
        'src/gene_tpm.R'

# 06 calculate cutoffs
rule calculate_cutoffs:
    input:
        bamfiles = expand(
            ('output/030_star-pass2/'
             '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
            stage=['UNM', 'PUNM', 'BCP', 'TCP'],
            plant=['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']),
        bg_counts = expand(
            ('output/050_calculate-background/'
             '{stage}_{plant}.csv'),
            stage=['UNM', 'PUNM', 'BCP', 'TCP'],
            plant=['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']),
    params:
        star_dir = 'output/030_star-pass2'
    output:
        intergenic_tpm = 'output/060_cutoffs/intergenic_tpm.csv'
    threads:
        1
    log:
        log = 'output/logs/060_cutoffs/intergenic_tpm.log'
    singularity:
        r_container
    script:
        'src/intergenic_tpm.R'


# 05. count reads per region
rule intergenic_reads:
    input:
        bam = ('output/030_star-pass2/'
               '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
        regions = 'output/040_shuffle/shuffled.gff3'
    output:
        counts = ('output/050_calculate-background/'
                  '{stage}_{plant}.csv')
    log:
        log = ('output/logs/050_calculate-background/'
               '{stage}-{plant}.log')
    threads:
        1
    singularity:
        r_container
    script:
        'src/count_reads_per_region.R'

# 4. prepare shuffle reference
rule shuffle:
    input:
        gtf = ('output/010_ref/'
               'Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf'),
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff',
        seqlengths = 'output/040_shuffle/seqlen.txt'
    output:
        shuffled = 'output/040_shuffle/shuffled.gff3'
    threads:
        1
    log:
        log = 'output/logs/040_shuffle/shuffle_gtf.log'
    singularity:
        r_container
    script:
        'src/shuffle_gtf.R'

rule calculate_seqlens:
    input:
        fasta = 'data/ref/TAIR10_Chr.all.fasta'
    output:
        seqlen = 'output/040_shuffle/seqlen.txt'
    threads:
        1
    singularity:
        bio_container
    script:
        'src/get_seqlength.py'

# 3. map
rule star_second_pass:
    input:
        r1 = 'output/020_trim-clip/{stage}_{plant}.fq.gz',
        star_reference = 'output/010_ref/star_reference/Genome',
        junctions = expand(
            'output/030_star-pass1/{stage}_{plant}.SJ.out.tab',
            stage=['UNM', 'PUNM', 'BCP', 'TCP'],
            plant=['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']),
    output:
        bam = ('output/030_star-pass2/'
               '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
        counts = 'output/030_star-pass2/{stage}_{plant}.ReadsPerGene.out.tab'
    threads:
        32
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/030_star-pass2/{stage}_{plant}.'
    log:
        'output/logs/030_STAR-pass2_{stage}_{plant}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outBAMcompression 10 '
        '--outReadsUnmapped Fastx '
        '--quantMode GeneCounts '
        '--readFilesIn {input.r1} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_first_pass:
    input:
        r1 = 'output/020_trim-clip/{stage}_{plant}.fq.gz',
        star_reference = 'output/010_ref/star_reference/Genome'
    output:
        sjdb = 'output/030_star-pass1/{stage}_{plant}.SJ.out.tab'
    threads:
        32
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/030_star-pass1/{stage}_{plant}.'
    log:
        'output/logs/030_STAR-pass1_{stage}-{plant}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '
        '--readFilesIn {input.r1} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

# 02 trim and clip with bbduk
rule trim_clip:
    input:
        unpack(find_input_files),
        adaptors = bbduk_adaptors,
        contaminants = bbduk_contaminants
    output:
        r1 = 'output/020_trim-clip/{stage}_{plant}.fq.gz'
    log:
        trim_log = 'output/logs/020_trim-clip/{stage}_{plant}_trim.log',
        trim_stats = 'output/020_trim-clip/{stage}_{plant}_trim-stats.txt',
        filter_log = 'output/logs/020_trim-clip/{stage}_{plant}_filter.log',
        filter_stats = 'output/020_trim-clip/{stage}_{plant}_filter-stats.txt'
    threads:
        1
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in={input.r1} '
        'out=stdout.fastq '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'ref={input.adaptors} '
        'stats={log.trim_stats} '
        'statscolumns=5 '
        '2> {log.trim_log} '
        '| bbduk.sh '
        'threads={threads} '
        '-Xmx100g '
        'in=stdin.fastq '
        'interleaved=f '
        'out={output.r1} '
        'ref={input.contaminants} '
        'k=31 hdist=1 stats={log.filter_stats} '
        'ziplevel=9 '
        'minlength=50 '
        '2> {log.filter_log}'

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
        32
    log:
        'output/logs/010_ref/star_reference.log'
    singularity:
        star_container
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
    singularity:
        r_container
    script:
        'src/preprocess_gtf.R'


