# Uncompress assembly


rule get_assembly:
    input:
        fa=lambda wc: MANIFEST[wc.sm]["asm"],
    output:
        fa=temp(join(OUTPUT_DIR, "mrsfast", "{sm}.fa")),
        fai=temp(join(OUTPUT_DIR, "mrsfast", "{sm}.fa.fai")),
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "jellyfish", "{sm}_get_assembly.log"),
    shell:
        """
        if [[ {input.fa} == *.gz ]]; then
            zcat {input.fa} > {output.fa} 2> {log}
        else
            ln -s {input.fa} {output.fa} 2> {log}
        fi
        samtools faidx {output.fa} 2> {log}
        """


# Get the kmers from the assembly
rule jellyfish_count:
    input:
        asm=rules.get_assembly.output.fa,
    output:
        counts=join(OUTPUT_DIR, "jellyfish", "{sm}_sunk.counts"),
    resources:
        mem=config["define_sunks"]["mem_jellyfish"],
    threads: config["define_sunks"]["threads_jellyfish"]
    conda:
        "../envs/viz.yaml"
    params:
        hash_size=10_000_000,
        upper_count=1,
        counter_len_bits=1,
        sunk_len=config["define_sunks"]["sunk_len"],
    log:
        join(LOG_DIR, "jellyfish", "{sm}_jellyfish_count.log"),
    benchmark:
        join(BMK_DIR, "jellyfish", "{sm}_jellyfish_count.tsv")
    shell:
        """
        jellyfish count \
        -m {params.sunk_len} \
        -s {params.hash_size} \
        -t {threads} \
        -C \
        -c {params.counter_len_bits} \
        -U {params.upper_count} {input.asm} \
        -o {output.counts} 2> {log}
        """


# Filter to kmers seen once in the assembly i.e. SUNKs
rule define_SUNKs:
    input:
        counts=rules.jellyfish_count.output.counts,
    output:
        db=join(OUTPUT_DIR, "jellyfish", "{sm}_sunk.db"),
        fa=join(OUTPUT_DIR, "jellyfish", "{sm}_sunk.fa"),
    resources:
        mem=config["define_sunks"]["mem_jellyfish"],
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "jellyfish", "{sm}_define_sunks.log"),
    benchmark:
        join(BMK_DIR, "jellyfish", "{sm}_define_sunks.tsv")
    threads: config["define_sunks"]["threads_jellyfish"]
    shell:
        """
        jellyfish dump -c -t {input.counts} | awk '{{print $1}}' > {output.db}
        awk '{{ print ">"$0"\\n"$0 }}' {output.db} > {output.fa}
        """


# Create index of assembly for mrsfast mapping
rule mrsfast_index:
    input:
        asm=rules.get_assembly.output.fa,
    output:
        index=temp(join(OUTPUT_DIR, "mrsfast", "{sm}.fa.index")),
    resources:
        mem=config["define_sunks"]["mem_mrsfast"],
    params:
        # Index with this window size.
        window_size=14,
    threads: config["define_sunks"]["threads_mrsfast"]
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "mrsfast", "{sm}_index.log"),
    benchmark:
        join(BMK_DIR, "mrsfast", "{sm}_index.tsv")
    shell:
        """
        mrsfast --ws {params.window_size} --index {input.asm} &> {log}
        """


# Map SUNKs back to assembly
rule mrsfast_search:
    input:
        asm=rules.get_assembly.output.fa,
        index=rules.mrsfast_index.output.index,
        db=rules.define_SUNKs.output.fa,
    output:
        sam=temp(join(OUTPUT_DIR, "mrsfast", "{sm}_sunk.sam")),
        bam=temp(join(OUTPUT_DIR, "mrsfast", "{sm}_sunk.bam")),
    resources:
        mem=config["define_sunks"]["mem_mrsfast"],
    params:
        # Need perfect match.
        err_thr=0,
    threads: config["define_sunks"]["threads_mrsfast"]
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "mrsfast", "{sm}_search.log"),
    benchmark:
        join(BMK_DIR, "mrsfast", "{sm}_search.tsv")
    shell:
        """
        mrsfast \
        --search {input.asm} \
        --threads {threads} \
        --mem {resources.mem} \
        --seq {input.db} \
        -o {output.sam} \
        --disable-nohits \
        -e {params.err_thr} &> {log}
        samtools sort -@ {threads} {output.sam} -o {output.bam} 2> {log}
        """


# Convert mrsfast bam file to bed file
rule bed_convert:
    input:
        bam=rules.mrsfast_search.output.bam,
    output:
        bed=temp(join(OUTPUT_DIR, "mrsfast", "{sm}_sunk.bed")),
        bedmerge=temp(join(OUTPUT_DIR, "mrsfast", "{sm}_sunk_merge.bed")),
        locs=join(OUTPUT_DIR, "mrsfast", "{sm}_sunk.loc"),
    resources:
        mem=16,
    conda:
        "../envs/viz.yaml"
    log:
        join(LOG_DIR, "mrsfast", "{sm}_bed_convert.log"),
    benchmark:
        join(BMK_DIR, "mrsfast", "{sm}_bed_convert.tsv")
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output.bed} 2> {log}
        bedtools merge -i {output.bed} > {output.bedmerge} 2> {log}
        {{ bedtools intersect -a {output.bed} -b {output.bedmerge} -wo | \
        cut -f 1,2,4,8 ;}} > {output.locs} 2> {log}
        """


rule define_sunks_all:
    input:
        expand(rules.bed_convert.output, sm=SAMPLES),
