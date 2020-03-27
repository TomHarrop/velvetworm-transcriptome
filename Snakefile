#!/usr/bin/env python3

from pathlib import Path
import multiprocessing
import pandas


#############
# FUNCTIONS #
#############

def get_lanes_for_indiv(wildcards):
    my_sd = raw_read_df.loc[(wildcards.indiv)]
    my_lanes = sorted(set(x[0] for x in my_sd.index.values))
    return snakemake.io.expand(
        'output/000_tmp/{indiv}_{lane}.repair.fastq',
        indiv=wildcards.indiv,
        lane=my_lanes)


def pick_trinity_input(wildcards):
    if wildcards.run == 'merged':
        return {
            'r1': expand(('output/000_tmp/'
                          '{indiv}.r1_joined_with_merged.fastq.gz'),
                         indiv=all_indivs),
            'r2': expand('output/020_merged/{indiv}_R2.fastq.gz',
                         indiv=all_indivs)}
    elif wildcards.run == 'raw':
        return {
            'r1': expand('output/010_reads/{indiv}_R1.fastq.gz',
                         indiv=all_indivs),
            'r2': expand('output/010_reads/{indiv}_R2.fastq.gz',
                         indiv=all_indivs)}
    else:
        raise ValueError(f'wtf run {wildcards.run}')


def resolve_input_fastq(wildcards):
    return raw_read_df.loc[(wildcards.indiv,
                            wildcards.lane,
                            wildcards.r),
                           'path']

###########
# GLOBALS #
###########

bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
trinity = 'shub://TomHarrop/assemblers:trinity_2.9.1'

raw_read_csv = 'data/raw_read_paths.csv'

########
# MAIN #
########

raw_read_df = pandas.read_csv(raw_read_csv,
                              index_col=['indiv', 'lane', 'read'])

all_indivs = sorted(set(x[0] for x in raw_read_df.index.values))
all_lanes = sorted(set(x[1] for x in raw_read_df.index.values))

#########
# RULES #
#########

rule target:
    input:
        expand('output/030_trinity/trinity_{run}/Trinity.fasta',
               run=['raw', 'merged'])

# trinity
rule trinity:
    input:
        unpack(pick_trinity_input)
    output:
        'output/030_trinity/trinity_{run}/Trinity.fasta',
        'output/030_trinity/trinity_{run}/Trinity.fasta.gene_trans_map'
    params:
        r1_line = lambda wildcards, input:
            ','.join(input.r1),
        r2_line = lambda wildcards, input:
            ','.join(input.r2),
        outdir = 'output/030_trinity/trinity_{run}'
    log:
        'output/logs/trinity.{run}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        trinity
    shell:
        'Trinity '
        # '--FORCE '
        '--seqType fq '
        '--max_memory 800G '
        '--left {params.r1_line} '
        '--right {params.r2_line} '
        '--SS_lib_type RF '
        '--CPU {threads} '
        '--output {params.outdir} '
        '&> {log}'

# merge the input reads, try with and without
rule join_merged_with_r1:
    input:
        r1 = 'output/020_merged/{indiv}_R1.fastq.gz',
        merged = 'output/020_merged/{indiv}_merged.fastq.gz'
    output:
        temp('output/000_tmp/{indiv}.r1_joined_with_merged.fastq.gz')
    singularity:
        bbduk
    shell:
        'cat {input.r1} {input.merged} > {output}'

rule merge:
    input:
        r1 = 'output/010_reads/{indiv}_R1.fastq.gz',
        r2 = 'output/010_reads/{indiv}_R2.fastq.gz'
    output:
        merged = 'output/020_merged/{indiv}_merged.fastq.gz',
        r1 = 'output/020_merged/{indiv}_R1.fastq.gz',
        r2 = 'output/020_merged/{indiv}_R2.fastq.gz',
        ihist = 'output/020_merged/{indiv}.ihist.txt'
    params:
        adaptors = '/adapters.fa'
    log:
        'output/logs/merge.{indiv}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbduk
    shell:
        'bbmerge.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.merged} '
        'outu1={output.r1} '
        'outu2={output.r2} '
        'ihist={output.ihist} '
        'verystrict=t '
        'adapters={params.adaptors} '
        '2>{log}'

# trinity doesn't do interleaved
rule split:
    input:
        'output/000_tmp/{indiv}.trim.fastq'
    output:
        r1 = 'output/010_reads/{indiv}_R1.fastq.gz',
        r2 = 'output/010_reads/{indiv}_R2.fastq.gz'
    log:
        'output/logs/split.{indiv}.log'
    singularity:
        bbduk
    shell:
        'reformat.sh '
        'in={input} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        'zl=9 '
        '2> {log}'

rule trim:
    input:
        'output/000_tmp/{indiv}.decon.fastq'
    output:
        pipe('output/000_tmp/{indiv}.trim.fastq')
    params:
        trim = '/adapters.fa',
        int_line = 'int=t '
    log:
        log = 'output/logs/trim.{indiv}.log',
        stats = 'output/logs/trim.{indiv}.stats'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'zl=9 '
        'in={input} '
        '{params.int_line} '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '


rule decon:
    input:
        'output/000_tmp/{indiv}.fastq'
    output:
        pipe('output/000_tmp/{indiv}.decon.fastq')
    params:
        filter = '/phix174_ill.ref.fa.gz',
        int_line = 'int=t '
    log:
        log = 'output/logs/decon.{indiv}.log',
        stats = 'output/logs/decon.{indiv}.stats'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'in={input} '
        '{params.int_line} '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '


rule combine:
    input:
        get_lanes_for_indiv
    output:
        pipe('output/000_tmp/{indiv}.fastq')
    singularity:
        bbduk
    shell:
        'cat {input} > {output}'

rule repair:
    input:
        r1 = 'output/000_tmp/{indiv}_{lane}_R1.barcode.fastq',
        r2 = 'output/000_tmp/{indiv}_{lane}_R2.barcode.fastq'
    output:
        pipe('output/000_tmp/{indiv}_{lane}.repair.fastq')
    log:
        'output/logs/repair.{indiv}_{lane}.log'
    singularity:
        bbduk
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'repair=t '
        '>> {output} '
        '2> {log}'


rule check_barcodes:
    input:
        resolve_input_fastq
    output:
        pipe('output/000_tmp/{indiv}_{lane}_{r}.barcode.fastq')
    params:
        bc = lambda wildcards:
            raw_read_df.loc[(wildcards.indiv,
                             wildcards.lane,
                             wildcards.r)].bc
    log:
        'output/logs/check_barcodes.{indiv}_{lane}_{r}.log'
    singularity:
        bbduk
    shell:
        'reformat.sh '
        'in={input} '
        'int=f '
        'out=stdout.fastq '
        'barcodefilter=t '
        'barcodes=\'{params.bc}\' '
        '>> {output} '
        '2> {log}'
