include: "gatherSplits.smk"


scattergather:
    split=config.get("nchunks", 10),


rule split_ONT:  # accept FOFN
    input:
        reads=lambda wc: MANIFEST[str(wc.sm)]["ont"],
    output:
        reads=temp(
            scatter.split(join(OUTPUT_DIR, "tag_ont", "{{sm}}_{scatteritem}.fq.gz"))
        ),
    resources:
        mem=2,
    threads: 8
    conda:
        "../envs/viz.yaml"
    log:
        scatter.split(join(LOG_DIR, "tag_ont", "{{sm}}_{scatteritem}_split_ONT.log")),
    shell:
        """
        cat {input.reads} | seqtk seq -F '#' | rustybam fastq-split {output.reads}
        """


# Determine sunk position on ONT reads.
# TODO: Can definitely rewrite to be faster. Single threaded and pre-compiled binary.
# Also why not use mrsfast here?
rule SUNK_annot:
    input:
        bin="workflow/scripts/kmerpos_annot3",
        locs=rules.bed_convert.output.locs,
        db=rules.define_SUNKs.output.db,
        ONT=join(OUTPUT_DIR, "tag_ont", "{sm}_{scatteritem}.fq.gz"),
    output:
        temp(join(OUTPUT_DIR, "tag_ont", "{sm}_{scatteritem}.sunkpos")),
    resources:
        mem=8,
    threads: 2
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "tag_ont", "{sm}_{scatteritem}_SUNK_annot.log"),
    benchmark:
        join(BMK_DIR, "tag_ont", "{sm}_{scatteritem}_SUNK_annot.tsv")
    shell:
        """
        {input.bin} {input.ONT} {input.db} {input.locs} {output} 2> {log}
        """


rule combine_ont_nofilt:
    input:
        gather_ONT_pos=gather.split(
            join(OUTPUT_DIR, "tag_ont", "{{sm}}_{scatteritem}.sunkpos")
        ),
    output:
        ONT_pos=join(OUTPUT_DIR, "tag_ont", "{sm}_detailed.sunkpos"),
    resources:
        mem=8,
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "tag_ont", "{sm}_combine_ont_nofilt.log"),
    shell:
        """
        cat {input.gather_ONT_pos} > {output.ONT_pos}
        """


# TODO: Why is this necessary? Just use samtools faidx.
rule read_lengths:
    input:
        ONT=join(OUTPUT_DIR, "tag_ont", "{sm}_{scatteritem}.fq.gz"),
    output:
        temp(join(OUTPUT_DIR, "tag_ont", "{sm}_{scatteritem}.rlen")),
    resources:
        mem=8,
        load=100,
    threads: 2
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "tag_ont", "{sm}_{scatteritem}_read_lengths.log"),
    benchmark:
        join(BMK_DIR, "tag_ont", "{sm}_{scatteritem}_read_lengths.tsv")
    shell:
        """
        workflow/scripts/rlen {input.ONT} {output} 2> {log}
        """


rule diag_filter_step:
    input:
        bin="workflow/scripts/diag_filter_v3",
        ONT_pos=rules.SUNK_annot.output,
        fai=rules.get_assembly.output.fai,
    output:
        ONT_pos_diag=temp(
            join(OUTPUT_DIR, "tag_ont", "{sm}_{scatteritem}_diag.sunkpos")
        ),
    resources:
        mem=8,
        load=100,
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "tag_ont", "{sm}_{scatteritem}_diag_filter_step.log"),
    benchmark:
        join(BMK_DIR, "tag_ont", "{sm}_{scatteritem}_diag_filter_step.tsv")
    shell:
        """
        {input.bin} {input.ONT_pos} {input.fai} > {output.ONT_pos_diag} 2> {log}
        """


rule diag_filter_final:
    input:
        bin="workflow/scripts/diag_filter_step2",
        ONT_pos=rules.SUNK_annot.output,
        ONT_pos_diag=rules.diag_filter_step.output.ONT_pos_diag,
    output:
        ONT_pos_diag_final=temp(
            join(OUTPUT_DIR, "tag_ont", "{sm}_{scatteritem}_diag2.sunkpos")
        ),
    resources:
        mem=8,
        load=100,
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "tag_ont", "{sm}_{scatteritem}_diag_filter_final.log"),
    benchmark:
        join(BMK_DIR, "tag_ont", "{sm}_{scatteritem}_diag_filter_final.tsv")
    shell:
        """
        {input.bin} {input.ONT_pos} {input.ONT_pos_diag} > {output.ONT_pos_diag_final} 2> {log}
        """


rule combine_ont:
    input:
        gather_ONT_pos=gather.split(
            join(OUTPUT_DIR, "tag_ont", "{{sm}}_{scatteritem}_diag2.sunkpos")
        ),
        gather_ONT_len=gather.split(
            join(OUTPUT_DIR, "tag_ont", "{{sm}}_{scatteritem}.rlen")
        ),
    output:
        ONT_pos=join(OUTPUT_DIR, "tag_ont", "{sm}.sunkpos"),
        ONT_len=join(OUTPUT_DIR, "tag_ont", "{sm}.rlen"),
    resources:
        mem=8,
        load=100,
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "tag_ont", "{sm}_combine_ont.log"),
    shell:
        """
        cat {input.gather_ONT_pos} > {output.ONT_pos}
        cat {input.gather_ONT_len} > {output.ONT_len}
        """


rule tag_ont_all:
    input:
        expand(rules.combine_ont_nofilt.output, sm=SAMPLES),
        expand(rules.combine_ont.output, sm=SAMPLES),
