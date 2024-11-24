

rule bad_sunks:
    input:
        script="workflow/scripts/badsunks_AR.py",
        fai=rules.get_assembly.output.fai,
        sunkpos=rules.combine_ont.output.ONT_pos,
    output:
        badsunks=join(OUTPUT_DIR, "final_out", "{sm}_bad_sunks.txt"),
    resources:
        mem=10,
        load=100,
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "final_out", "{sm}_bad_sunks.log"),
    benchmark:
        join(BMK_DIR, "final_out", "{sm}_bad_sunks.tsv")
    shell:
        """
        python {input.script} {input.fai} {input.sunkpos} {output.badsunks}
        """


checkpoint split_sunkpos:
    input:
        script="workflow/scripts/split_locs.py",
        ONT_pos=rules.combine_ont.output.ONT_pos,
        kmer_loc=rules.bed_convert.output.locs,
    output:
        temp(directory(join(OUTPUT_DIR, "final_out", "{sm}_breaks"))),
    resources:
        mem=10,
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "final_out", "{sm}_split_sunkpos.log"),
    shell:
        """
        python {input.script} --ont_pos {input.ONT_pos} --kmer_loc {input.kmer_loc} --outdir {output} 2> {log}
        """


rule process_by_contig:
    input:
        script="workflow/scripts/process-by-contig_lowmem_AR.py",
        sunk_pos_contig=join(OUTPUT_DIR, "final_out", "{sm}_breaks", "{contig}.sunkpos"),
        locs_contig=join(OUTPUT_DIR, "final_out", "{sm}_breaks", "{contig}.loc"),
        rlen=rules.combine_ont.output.ONT_len,
        bad_sunks=rules.bad_sunks.output.badsunks,
    output:
        outputdf=join(OUTPUT_DIR, "final_out", "{sm}_inter_outs", "{contig}.tsv"),
        bed=join(OUTPUT_DIR, "final_out", "{sm}_bed_files", "{contig}.bed"),
    resources:
        mem=16,
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "final_out", "{sm}_{contig}_process_by_contig.log"),
    benchmark:
        join(BMK_DIR, "final_out", "{sm}_{contig}_process_by_contig.tsv")
    shell:
        """
        python {input.script} \
        {input.locs_contig} \
        {input.sunk_pos_contig} \
        {input.rlen} \
        {input.bad_sunks} \
        {output.outputdf} \
        {output.bed} 2> {log}
        """


rule gather_process_by_contig:
    input:
        bed=gatherAsmBeds,
    output:
        allbeds=join(OUTPUT_DIR, "final_out", "{sm}.valid.bed"),
    resources:
        mem=8,
        load=100,
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "final_out", "{sm}_gather_process_by_contig.log"),
    shell:
        """
        cat {input.bed} > {output.allbeds}
        """


# TODO: This isn't great. Why not just past the merged bed file?
rule get_gaps:
    input:
        script="workflow/scripts/get_gaps.py",
        fai=rules.get_assembly.output.fai,
        allbeds=rules.gather_process_by_contig.output.allbeds,
    output:
        gaps=join(OUTPUT_DIR, "final_out", "{sm}.gaps.bed"),
        nodata=join(OUTPUT_DIR, "final_out", "{sm}.nodata.bed"),
    resources:
        mem=10,
    params:
        # Expects {input_dir}/{contig}.bed
        input_dir=lambda wc: f"results/{wc.sm}/bed_files/",
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "final_out", "{sm}_get_gaps.log"),
    shell:
        """
        python {input.script} \
        {input.fai} \
        {params.input_dir} \
        {output.gaps} {output.nodata} 2> {log}
        """


rule slop_gaps:
    input:
        gaps=rules.get_gaps.output.gaps,
        fai=rules.get_assembly.output.fai,
    output:
        gaps_slop=join(OUTPUT_DIR, "final_out", "{sm}.gaps.slop.bed"),
    resources:
        mem=10,
    params:
        bp=200_000,
    threads: 1
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "final_out", "{sm}_slop_gaps.log"),
    shell:
        """
        bedtools slop -i {input.gaps} -g {input.fai} -b {params.bp} > {output.gaps_slop}
        """


rule validate_sunks_all:
    input:
        expand(rules.bad_sunks.output, sm=SAMPLES),
        expand(rules.gather_process_by_contig.output, sm=SAMPLES),
        expand(rules.get_gaps.output, sm=SAMPLES),
        expand(rules.slop_gaps.output, sm=SAMPLES),
