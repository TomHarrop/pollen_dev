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
r_container = ('shub://TomHarrop/singularity-containers:r_3.5.1'
               '@6fa6a4a5a2b6da669923db2a3b8a0bb3f876003c')
star_container = ('shub://TomHarrop/singularity-containers:star_2.6.0c'
                  '@eaa90a258fdb26b6b0ce7d07246ffe2c')
bbduk_container = ('shub://TomHarrop/singularity-containers:bbmap_38.00'
                   '@a773baa8cc025cc5b5cbee20e507fef7')
biop_container = ('shub://TomHarrop/singularity-containers:biopython_1.72'
                  '@eb4213531ecbfb44bc04028e2c3b1559')
bioc_container = ('shub://TomHarrop/singularity-containers:bioconductor_3.7'
                  '@2785c89cc4ef1cbb08c06dde9ecf9544')

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
        'output/090_deseq/dds.Rds',
        'output/070_tpm/tpm_summary.csv',
        'output/090_deseq/pca.csv',
        'output/090_deseq/wald_stage.csv',
        'output/050_calculate-background/featurecounts.csv',
        'output/100_venn-diagrams/array_comparison.csv'

# 10 set analysis
rule venn_diagram:
    input:
        array = 'data/gb-2004-5-11-r85-s1.xls',
        calls = 'output/080_filter-background/gene_calls.csv'
    output:
        venn_diagram = 'output/100_venn-diagrams/venn_diagram.pdf',
        array_comparison = 'output/100_venn-diagrams/array_comparison.csv'
    params:
        outdir = 'output/100_venn-diagrams'
    log:
        'output/logs/100_venn/venn_diagram.log'
    threads:
        1
    singularity:
        r_container
    script:
        'src/venn_diagram.R'


# 09 DESeq analysis
rule deseq_wald_tests:
    input:
        dds = 'output/090_deseq/dds.Rds'
    output:
        group_test = 'output/090_deseq/wald_LBCP-LTCP_vs_RUNM-PUNM.csv',
        stage_tests = 'output/090_deseq/wald_stage.csv',
        de_matrix = 'output/090_deseq/number_of_de_genes.csv'
    params:
        alpha = 0.1,
        lfc_threshold = 0.5849625       # log(1.5, 2)
    threads:
        20
    log:
        log = 'output/logs/090_deseq/deseq_wald_tests.log'
    singularity:
        bioc_container
    script:
        'src/deseq_wald_tests.R'

rule deseq_qc:
    input:
        dds = 'output/090_deseq/dds.Rds'
    output:
        pca_plot = 'output/090_deseq/pca_plot.pdf',
        pca_dt = 'output/090_deseq/pca.csv',
        distance_heatmap = 'output/090_deseq/distance_heatmap.pdf'
    threads:
        20
    log:
        log = 'output/logs/090_deseq/deseq_qc.log'
    singularity:
        bioc_container
    script:
        'src/deseq_qc.R'

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
        bioc_container
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
rule summarise_tpm:
    input:
        tpm = 'output/080_filter-background/gene_calls.csv'
    output:
        tpm_wide = 'output/070_tpm/tpm_wide.csv',
        tpm_summary = 'output/070_tpm/tpm_summary.csv',
        tpm_summary_wide = 'output/070_tpm/tpm_summary_wide.csv'
    threads:
        1
    log:
        log = 'output/logs/070_tpm/summarise_tpm.log'
    singularity:
        r_container
    script:
        'src/summarise_tpm.R'

rule calculate_tpm:
    input:
        count_files = expand(
            ('output/030_star-pass2/{stage}_{plant}.ReadsPerGene.out.tab'),
            stage=['RUNM', 'PUNM', 'LBCP', 'LTCP'],
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
        bioc_container
    script:
        'src/gene_tpm.R'

# 06 calculate cutoffs
rule calculate_cutoffs:
    input:
        bamfiles = expand(
            ('output/030_star-pass2/'
             '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
            stage=['RUNM', 'PUNM', 'LBCP', 'LTCP'],
            plant=['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']),
        bg_counts = expand(
            ('output/050_calculate-background/'
             '{stage}_{plant}.csv'),
            stage=['RUNM', 'PUNM', 'LBCP', 'LTCP'],
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
rule feature_counts:
    input:
        bam_files = expand(
            ('output/030_star-pass2/'
             '{stage}_{plant}.Aligned.sortedByCoord.out.bam'),
            stage=['RUNM', 'PUNM', 'LBCP', 'LTCP'],
            plant=['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']),
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff'
    output:
        feature_counts = 'output/050_calculate-background/featurecounts.csv'
    log:
        log = 'output/logs/050_calculate-background/feature_counts.log'
    threads:
        50
    singularity:
        bioc_container
    script:
        'src/count_reads_per_feature.R'

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
        bioc_container
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
        bioc_container
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
        biop_container
    script:
        'src/get_seqlength.py'

# 3. map
rule star_second_pass:
    input:
        r1 = 'output/020_trim-clip/{stage}_{plant}.fq.gz',
        star_reference = 'output/010_ref/star_reference/Genome',
        junctions = expand(
            'output/030_star-pass1/{stage}_{plant}.SJ.out.tab',
            stage=['RUNM', 'PUNM', 'LBCP', 'LTCP'],
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

# parse annotations from araport GFF
rule make_annotation_table:
    input:
        gff = 'data/ref/Araport11_GFF3_genes_transposons.201606.gff'
    output:
        annotation = 'output/010_ref/araport_annotation.csv'
    threads:
        1
    log:
        'output/logs/010_ref/make_annotation_table.log'
    singularity:
        bioc_container
    script:
        'src/make_annotation_table.R'


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
        bioc_container
    script:
        'src/preprocess_gtf.R'


