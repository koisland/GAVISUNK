def gatherAsmBeds(wc):
    breaks_dir = checkpoints.split_sunkpos.get(**wc).output[0]
    contigs = glob_wildcards(join(breaks_dir, "{contig}.sunkpos")).contig
    return expand(
        rules.process_by_contig.output.bed,
        sm=wc.sm,
        contig=contigs,
    )


def getSunkLocs(wc):
    breaks_dir = checkpoints.split_sunkpos.get(**wc).output[0]
    contigs = glob_wildcards(join(breaks_dir, "{contig}.loc")).contig
    return expand(
        "results/{sample}/breaks/{contigs}_{hap}.loc",
        sm=wc.sm,
        hap=wildcards.hap,
        contigs=contigs,
    )


def getInterOut(wildcards):
    CONTIGS = glob_wildcards(
        "results/{sample}/inter_outs/{contigs}_{hap}.tsv".format(
            sample=wildcards.sample, hap=wildcards.hap, contigs="{contigs}"
        )
    ).contigs
    # print("results/{sample}/inter_outs/{contigs}_{hap}.tsv")
    return expand(
        "results/{sample}/inter_outs/{contigs}_{hap}.tsv",
        sample=wildcards.sample,
        hap=wildcards.hap,
        contigs=CONTIGS,
    )


def getSunkPos(wildcards):
    CONTIGS = glob_wildcards(
        "results/{sample}/breaks/{contigs}_{hap}.sunkpos".format(
            sample=wildcards.sample, hap=wildcards.hap, contigs="{contigs}"
        )
    ).contigs
    return expand(
        "results/{sample}/breaks/{contigs}_{hap}.sunkpos",
        sample=wildcards.sample,
        hap=wildcards.hap,
        contigs=CONTIGS,
    )


def vizInputs(wildcards):
    # print(wildcards.runmode)
    if wildcards.runmode == "user_bed":
        bed = manifest_df.at[wildcards.sample, f"{wildcards.hap}_bed"]
    elif wildcards.runmode == "gaps":
        bed = rules.slop_gaps.output.gaps_slop
    input_dict = {
        "bed": bed,
        "rlen": rules.combine_ont.output.ONT_len,
        "interout": rules.confirm_out.output.flag,
        "pos_locs": rules.split_sunkpos.output.flag,
    }
    if not pd.isnull(manifest_df.at[wildcards.sample, f"{wildcards.hap}_colortrack"]):
        input_dict["colorbed"] = manifest_df.at[
            wildcards.sample, f"{wildcards.hap}_colortrack"
        ]
    # print(input_dict)
    return input_dict


def vizInputsDetailed(wildcards):
    # print(wildcards.runmode)
    if wildcards.runmode == "user_bed":
        bed = manifest_df.at[wildcards.sample, f"{wildcards.hap}_bed"]
    elif wildcards.runmode == "gaps":
        bed = rules.slop_gaps.output.gaps_slop
    input_dict = {
        "bed": bed,
        "rlen": rules.combine_ont.output.ONT_len,
        "interout": rules.confirm_out.output.flag,
        "pos_locs": rules.combine_ont_nofilt.output.ONT_pos,
        "splits": rules.split_sunkpos.output.flag,
    }
    if not pd.isnull(manifest_df.at[wildcards.sample, f"{wildcards.hap}_colortrack"]):
        input_dict["colorbed"] = manifest_df.at[
            wildcards.sample, f"{wildcards.hap}_colortrack"
        ]
    return input_dict


def covprobInputs(wildcards):
    if wildcards.hap == "hap1":
        bed = rules.get_gaps.output.gaps1
    elif wildcards.hap == "hap2":
        bed = rules.get_gaps.output.gaps2
    else:
        raise ValueError(f"UNKNOWN HAP: {wildcards.hap}")
    input_dict = {
        "bed": bed,
        "locs": rules.bed_convert.output.locs,
        "rlen": rules.combine_ont.output.ONT_len,
    }
    return input_dict
