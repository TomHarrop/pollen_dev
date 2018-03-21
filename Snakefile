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

sample_key_file = 'data/sample_key.csv'
read_dir = 'data/fastq'
bbduk_adaptors = 'data/bbmap_resources/adapters.fa'
bbduk_contaminants = 'data/bbmap_resources/sequencing_artifacts.fa.gz'
star_reference_folder = 'output/010_ref/star_reference'

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
        expand(('output/030_star-pass2/'
                '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
               stage=['UNM', 'PUNM', 'BCP', 'TCP'],
               plant=['p1', 'p2', 'p3', 'p4']),
        expand(('output/050_calculate-background/{region}/'
                '{stage}_{plant}.csv'),
               stage=['UNM', 'PUNM', 'BCP', 'TCP'],
               plant=['p1', 'p2', 'p3', 'p4'],
               region=['genic', 'intergenic'])


# 05. count reads per region
rule intergenic_reads:
    input:
        bam = ('output/030_star-pass2/'
               '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
        regions = 'output/040_shuffle/shuffled.bed'
    output:
        counts = ('output/050_calculate-background/intergenic/'
                  '{stage}_{plant}.csv')
    threads:
        1
    script:
        'src/count_reads_per_region.R'

rule genic_reads:
    input:
        bam = ('output/030_star-pass2/'
               '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
        regions = 'output/040_shuffle/genes.bed'
    output:
        counts = ('output/050_calculate-background/genic/'
                  '{stage}_{plant}.csv')
    threads:
        1
    script:
        'src/count_reads_per_region.R'


# 4. prepare shuffle reference
rule shuffle:
    input:
        genic = 'output/040_shuffle/genes.bed',
        intergenic = 'output/040_shuffle/intergenic.bed',
        seqlen = 'output/040_shuffle/seqlen.txt'
    output:
        'output/040_shuffle/shuffled.bed'
    threads:
        1
    log:
        'output/logs/040_shuffle/shuffle.log'
    shell:
        'bedtools shuffle '
        '-incl {input.intergenic} '
        '-noOverlapping '
        '-seed 7 '
        '-i {input.genic} '
        '-g {input.seqlen} '
        '> {output} '
        '2> {log}'

rule split_gffs:
    input:
        gtf = ('output/010_ref/'
               'Araport11_GFF3_genes_transposons_nuc_norrna.201606.gtf'),
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff'
    output:
        genic = 'output/040_shuffle/genes.bed',
        intergenic = 'output/040_shuffle/intergenic.bed'
    threads:
        1
    log:
        log = 'output/logs/040_shuffle/split_gffs.log'
    script:
        'src/split_gtf.R'

rule calculate_seqlens:
    input:
        fasta = 'data/ref/TAIR10_Chr.all.fasta'
    output:
        seqlen = 'output/040_shuffle/seqlen.txt'
    threads:
        1
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
            plant=['p1', 'p2', 'p3', 'p4'])
    output:
        bam = ('output/030_star-pass2/'
               '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
        counts = 'output/030_star-pass2/{stage}_{plant}.ReadsPerGene.out.tab'
    threads:
        10
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/030_star-pass2/{stage}_{plant}.'
    log:
        'output/logs/030_STAR-pass2_{stage}_{plant}.log'
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
        10
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/030_star-pass1/{stage}_{plant}.'
    log:
        'output/logs/030_STAR-pass1_{stage}-{plant}.log'
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
        10
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


