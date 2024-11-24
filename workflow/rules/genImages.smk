include: "gatherSplits.smk"


rule viz_contigs_detailed:
    input:
        unpack(vizInputsDetailed),
    output:
        flag=touch("results/pngs/{runmode}/{sample}/{sample}_{hap}_detailed.done"),
    resources:
        mem=160,
        load=200,
    threads: 2
    conda:
        "../envs/viz.yaml"
    log:
        "logs/{runmode}_{sample}_{hap}_viz_contigs_detailed.log",
    script:
        "../scripts/viz_detailed.py"


rule viz_contigs:
    input:
        unpack(vizInputs),
    output:
        flag=touch("results/pngs/{runmode}/{sample}/{sample}_{hap}.done"),
    resources:
        mem=160,
        load=200,
    threads: 2
    conda:
        "../envs/viz.yaml"
    log:
        "logs/{runmode}_{sample}_{hap}_viz_contigs.log",
    script:
        "../scripts/viz.py"


rule covprob:
    input:
        unpack(covprobInputs),
        script="workflow/scripts/covprob.py",
        fai=rules.get_assembly.output.fai,
    output:
        tsv="results/{sample}/final_out/{hap}.gaps.covprob.tsv",
    resources:
        mem=160,
        load=200,
    threads: 2
    benchmark:
        "benchmarks/{sample}_{hap}_covprob.bench"
    conda:
        "../envs/viz.yaml"
    log:
        "logs/{sample}_{hap}_covprob.log",
    shell:
        """
        python {input.script}
        --fai {input.fai}
        --bed
        --locs
        --rlen
        --kmer_len
        --outfile
        """


# if BED == "":
#     bed_df = pd.DataFrame()
# else:
#     bed_df = pd.read_csv(
#         BED, sep="\t", header=None, names=["contig", "start", "stop", "gene"]
#     )
#     bed_df.set_index(["gene"], inplace=True, drop=True)


detailed_plot = config.get("plot_detailed", False)

if detailed_plot:
    pass
    # target = (
    #     expand(
    #         "results/pngs/gaps/{sample}/{sample}_{hap}_detailed.done",
    #         sample=manifest_df.index,
    #         hap=["hap1", "hap2"],
    #     ),
    # )
