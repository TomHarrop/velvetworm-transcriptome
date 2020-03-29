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
busco = 'docker://ezlabgva/busco:v4.0.4_cv1'
r = 'shub://TomHarrop/r-containers:r_3.6.2'
trinity = 'shub://TomHarrop/assemblers:trinity_2.9.1'
trinotate = 'shub://TomHarrop/trinotate_pipeline:v0.0.12'


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
        # expand('output/030_trinity/trinity_{run}/Trinity.fasta',
        #        run=['raw', 'merged'])
        expand('output/040_trinity-abundance/{run}/ExN50.pdf',
               run=['raw', 'merged']),
        expand('output/045_transcript-length/{run}/blastx.outfmt6.grouped',
               run=['raw', 'merged']),
        'output/099_busco/busco_plot.pdf'


# busco
rule plot_busco_results:
    input:
        busco_files = expand(
            'output/099_busco/{run}.{filter}/fixed_full_table.csv',
            run=['merged', 'raw'],
            filter=['length', 'expr'])
    output:
        plot = 'output/099_busco/busco_plot.pdf'
    log:
        'output/logs/plot_busco_results.log'
    singularity:
        r
    script:
        'src/plot_busco_results.R'

rule fix_busco_table:
    input:
        'output/099_busco/{run}.{filter}/busco/run_hymenoptera_odb10/full_table.tsv'
    output:
        'output/099_busco/{run}.{filter}/fixed_full_table.csv'
    log:
        'output/logs/fixed_full_table.{run}.{filter}.log'
    singularity:
        r
    script:
        'src/fix_busco_table.R'

rule busco:
    input:
        fasta = 'output/050_transcripts-by/{run}/{filter}.fasta',
        lineage = 'data/hymenoptera_odb10'
    output:
        ('output/099_busco/{run}.{filter}/'
         'busco/run_hymenoptera_odb10/'
         'full_table.tsv'),
    log:
        Path(('output/logs/'
              'busco.{run}.{filter}.log')).resolve()
    params:
        wd = 'output/099_busco/{run}.{filter}',
        name = 'busco',
        fasta = lambda wildcards, input: Path(input.fasta).resolve(),
        lineage = lambda wildcards, input:
            Path(input.lineage).resolve()
    threads:
        multiprocessing.cpu_count()
    singularity:
        busco
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.fasta} '
        '--out {params.name} '
        '--lineage_dataset {params.lineage} '
        '--cpu {threads} '
        '--augustus_species honeybee1 '
        '--mode transcriptome '
        '&> {log}'

# annotation
rule trinotate:
    input:
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
        blast_db = 'data/db/uniprot_sprot.pep',
        hmmer_db = 'data/db/Pfam-A.hmm',
        sqlite_db = 'data/db/Trinotate.sqlite',
    output:
        'output/070_trinotate/{run}/trinotate/trinotate_annotation_report.txt',
        'output/070_trinotate/{run}/blastx/blastx.outfmt6',
        'output/070_trinotate/{run}/trinotate/Trinotate.sqlite'
    params:
        outdir = 'output/070_trinotate/{run}'
    log:
        'output/logs/trinotate.{run}.log'
    threads:
        min(multiprocessing.cpu_count, 64)
    singularity:
        trinotate
    shell:
        'trinotate_pipeline '
        '--trinity_fasta {input.transcripts} '
        '--blast_db {input.blast_db} '
        '--hmmer_db {input.hmmer_db} '
        '--sqlite_db {input.sqlite_db} '
        '--outdir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


# analyse Trinity output
rule group_blast_hits:
    input:
        blastx = 'output/070_trinotate/{run}/blastx/blastx.outfmt6',
        db = 'data/db/uniprot_sprot.pep',
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta'
    output:
        'output/045_transcript-length/{run}/blastx.outfmt6.grouped',
        ('output/045_transcript-length/{run}/'
         'blastx.outfmt6.grouped.w_pct_hit_length.txt')
    params:
        wd = 'output/045_transcript-length/{run}',
        transcripts = lambda wildcards, input:
            Path(input.transcripts).resolve(),
        blastx = lambda wildcards, input:
            Path(input.blastx).resolve(),
        db = lambda wildcards, input:
            Path(input.db).resolve()
    log:
        Path('output/logs/group_blast_hits.{run}.log').resolve()
    singularity:
        trinity
    shell:
        'cd {params.wd} || exit 1 ; '
        'ln -s {params.blastx} blastx-group.outfmt6 ; '
        'blast_outfmt6_group_segments.pl '
        'blastx-group.outfmt6 '
        '{params.transcripts} '
        '{params.db} '
        '> blastx.outfmt6.grouped '
        '2> {log} ;'
        'blast_outfmt6_group_segments.tophit_coverage.pl '
        'blastx.outfmt6.grouped '
        '> blastx.outfmt6.grouped.w_pct_hit_length.txt '
        '2>> {log}'

rule add_blast_coverage:
    input:
        blastx = 'output/070_trinotate/{run}/blastx/blastx.outfmt6',
        db = 'data/db/uniprot_sprot.pep',
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta'
    output:
        ('output/045_transcript-length/{run}/'
         'blastx.outfmt6.w_pct_hit_length.txt'),
        ('output/045_transcript-length/{run}/'
         'blastx-coverage.outfmt6.w_pct_hit_length')
    params:
        wd = 'output/045_transcript-length/{run}',
        transcripts = lambda wildcards, input:
            Path(input.transcripts).resolve(),
        blastx = lambda wildcards, input:
            Path(input.blastx).resolve(),
        db = lambda wildcards, input:
            Path(input.db).resolve()
    log:
        Path('output/logs/add_blast_coverage.{run}.log').resolve()
    singularity:
        trinity
    shell:
        'cd {params.wd} || exit 1 ; '
        'ln -s {params.blastx} blastx-coverage.outfmt6 ; '
        'analyze_blastPlus_topHit_coverage.pl '
        'blastx-coverage.outfmt6 '
        '{params.transcripts} '
        '{params.db} '
        '> blastx.outfmt6.w_pct_hit_length.txt '
        '2> {log}'

rule plot_exn50:
    input:
        exn50 = 'output/040_trinity-abundance/{run}/ExN50.stats'
    output:
        plot = 'output/040_trinity-abundance/{run}/ExN50.pdf'
    log:
        'output/logs/plot_exn50.{run}.log'
    singularity:
        r
    script:
        'src/plot_exn50.R'

rule contig_exn50:
    input:
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
        expr = ('output/040_trinity-abundance/{run}/'
                'salmon.isoform.TMM.EXPR.matrix')
    output:
        inputs = ('output/040_trinity-abundance/{run}/'
                  'salmon.isoform.TMM.EXPR.matrix.E-inputs'),
        stats = ('output/040_trinity-abundance/{run}/'
                 'ExN50.stats')
    params:
        outdir = 'output/040_trinity-abundance/{run}',
        transcripts = lambda wildcards, input:
            Path(input.transcripts).resolve(),
        expr = lambda wildcards, input:
            Path(input.expr).resolve()
    log:
        Path('output/logs/abundance_to_matrix.{run}.log').resolve()
    singularity:
        trinity
    shell:
        'cd {params.outdir} || exit 1 ; '
        'contig_ExN50_statistic.pl '
        '{params.expr} '
        '{params.transcripts} '
        '> ExN50.stats '
        '2> {log}'


rule abundance_to_matrix:
    input:
        qf = expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
                    indiv=all_indivs),
        gtm = 'output/030_trinity/trinity_{run}/Trinity.fasta.gene_trans_map'
    output:
        'output/040_trinity-abundance/{run}/salmon.isoform.counts.matrix',
        'output/040_trinity-abundance/{run}/salmon.isoform.TMM.EXPR.matrix'
    params:
        outdir = 'output/040_trinity-abundance/{run}',
        qf = lambda wildcards, input:
            sorted(set(Path(x).resolve() for x in input.qf)),
        gtm = lambda wildcards, input:
            Path(input.gtm).resolve()
    log:
        Path('output/logs/abundance_to_matrix.{run}.log').resolve()
    singularity:
        trinity
    shell:
        'cd {params.outdir} || exit 1 ; '
        'abundance_estimates_to_matrix.pl '
        '--est_method salmon '
        '--gene_trans_map {params.gtm} '
        '--name_sample_by_basedir '
        '{params.qf} '
        '&> {log}'

rule trinity_abundance:
    input:
        'output/030_trinity/trinity_{run}/Trinity.fasta.salmon.idx',
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
        r1 = expand('output/010_reads/{indiv}_R1.fastq.gz',
                    indiv=all_indivs),
        r2 = expand('output/010_reads/{indiv}_R2.fastq.gz',
                    indiv=all_indivs)
    output:
        expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
               indiv=all_indivs)
    params:
        outdir = 'output/040_trinity-abundance/{run}',
        transcripts = lambda wildcards, input:
            Path(input.transcripts).resolve(),
        r1_line = lambda wildcards, input:
            ','.join(input.r1),
        r2_line = lambda wildcards, input:
            ','.join(input.r2),
    log:
        Path('output/logs/trinity_abundance.{run}.log').resolve()
    threads:
        multiprocessing.cpu_count()
    singularity:
        trinity
    shell:
        'cd {params.outdir} || exit 1 ; '
        'align_and_estimate_abundance.pl '
        '--transcripts {params.transcripts} '
        '--seqType fq '
        '--left {params.r1_line} '
        '--right {params.r2_line} '
        '--est_method salmon '
        '--SS_lib_type RF '
        '--thread_count {threads} '
        '--trinity_mode '
        '&> {log}'


rule trinity_abundance_prep:
    input:
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta'
    output:
        directory('output/030_trinity/trinity_{run}/Trinity.fasta.salmon.idx')
    log:
        'output/logs/trinity_abundance_prep.{run}.log'
    singularity:
        trinity
    shell:
        'align_and_estimate_abundance.pl '
        '--transcripts {input.transcripts} '
        '--est_method salmon '
        '--trinity_mode '
        '--prep_reference '
        '&> {log}'


# filter transcripts
rule filter_isoforms:
    input:
        transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
        names = 'output/050_transcripts-by/{run}/transcripts_by_{filter}.txt'
    output:
        'output/050_transcripts-by/{run}/{filter}.fasta'
    log:
        'output/logs/filter_isoforms.{run}.{filter}.log'
    singularity:
        bbduk
    shell:
        'filterbyname.sh '
        'in={input.transcripts} '
        'include=t '
        'names={input.names} '
        'out={output} '
        '2> {log}'

rule transcripts_by:
    input:
        qf = expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
                    indiv=all_indivs),
        expr = ('output/040_trinity-abundance/{run}/'
                'salmon.isoform.TMM.EXPR.matrix'),
        gtm = 'output/030_trinity/trinity_{run}/Trinity.fasta.gene_trans_map'
    output:
        ln = ('output/050_transcripts-by/{run}/'
              'transcripts_by_length.txt'),
        expr = ('output/050_transcripts-by/{run}/'
                'transcripts_by_expr.txt')
    log:
        'output/logs/transcripts_by.{run}.log'
    singularity:
        r
    script:
        'src/transcripts_by.R'


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
