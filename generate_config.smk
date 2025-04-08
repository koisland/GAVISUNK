
run_name = config["run_name"]

# WCS_READS = glob_wildcards("bed/{asm}/{sm}.{hap}.txt")
WCS_READS = glob_wildcards(f"bed/{run_name}/{{sm}}.{{hap}}.txt")
# ASM, SM, HAP = WCS_READS.asm, WCS_READS.sm, WCS_READS.hap 
ASM, SM, HAP = [run_name for _ in range(len(WCS_READS.sm))], WCS_READS.sm, WCS_READS.hap 
WCS = dict(asm=ASM, sm=SM, hap=HAP)


rule generate_subset_asm_fai:
    input:
        bed="bed/{asm}/{sm}.{hap}.txt",
        asm_fai="asm/{asm}/{sm}-asm-comb-dedup.fa.fai"
    output:
        fai="fai/{asm}/{sm}.{hap}.subset.fai"
    shell:
        """
        grep -f <(cut -f 1 {input.bed}) {input.asm_fai} > {output.fai}
        """

rule generate_manifest:
    input:
        cov="sm_cov.tsv",
        asm=expand("asm/{asm}/{sm}-asm-comb-dedup.fa", zip, **WCS),
        fai=expand(rules.generate_subset_asm_fai.output, zip, **WCS),
        reads=expand("reads/{asm}/{sm}.subset.fasta", zip, **WCS),
        bed=expand("bed/{asm}/{sm}.{hap}.txt", zip, **WCS),
    output:
        f"hgsvc_{run_name}.tsv"
    run:
        header=[
            "sample",
            "hap1_ONT",
            "hap2_ONT",
            "asm",
            "hap1_bed",
            "hap2_bed",
            "hap1_colortrack",
            "hap2_colortrack",
            "hap1_cov",
            "hap2_cov",
            "hap1_asm_fai",
            "hap2_asm_fai",
        ]
        with open(str(input.cov), "rt") as fh:
            sm_cov = {}
            for line in fh.readlines():
                sm, cov = line.strip().split("\t")
                sm_cov[sm] = float(cov)

        with open(str(output[0]), "wt") as fh:
            print("\t".join(header), file=fh)
            comb = set((asm, sm) for asm, sm in zip(ASM, SM))
            bedfiles = set(input.bed)
            for asm, sm in comb:
                cov = sm_cov[sm]
                hap_cov = cov // 2
                bed_hap_1 = f"bed/{asm}/{sm}.haplotype1.txt"
                bed_hap_2 = f"bed/{asm}/{sm}.haplotype2.txt"
                if bed_hap_1 in bedfiles and bed_hap_2 in bedfiles:
                    comb_bedfiles = [
                        f"bed/{asm}/{sm}.haplotype1.txt",
                        f"bed/{asm}/{sm}.haplotype2.txt",
                    ]
                    comb_fai = [
                        expand(rules.generate_subset_asm_fai.output, sm=sm, hap="haplotype1", asm=asm)[0],
                        expand(rules.generate_subset_asm_fai.output, sm=sm, hap="haplotype2", asm=asm)[0],
                    ]
                elif bed_hap_1 in bedfiles:
                    comb_bedfiles = [f"bed/{asm}/{sm}.haplotype1.txt"] * 2
                    comb_fai = [
                        expand(rules.generate_subset_asm_fai.output, sm=sm, hap="haplotype1", asm=asm)[0]
                    ] * 2
                elif bed_hap_2 in bedfiles:
                    comb_bedfiles = [f"bed/{asm}/{sm}.haplotype2.txt"] * 2
                    comb_fai = [
                        expand(rules.generate_subset_asm_fai.output, sm=sm, hap="haplotype2", asm=asm)[0]
                    ] * 2
                
                line = [
                    f"{sm}_{asm}",
                    f"reads/{asm}/{sm}.subset.fasta",
                    f"reads/{asm}/{sm}.subset.fasta",
                    f"asm/{asm}/{sm}-asm-comb-dedup.fa",
                    *comb_bedfiles,
                    "",
                    "",
                    str(hap_cov),
                    str(hap_cov),
                    *comb_fai,
                ]
                print("\t".join(line), file=fh)


rule generate_config:
    input:
        cfg="config/config.yaml",
        manifest=rules.generate_manifest.output
    output:
        cfg=f"config_hgsvc_{run_name}.yaml"
    params:
        sunk_len = config["sunk_len"]
    run:
        import yaml

        with (
            open(input.cfg, "rb") as fh,
            open(output.cfg, "wt") as out_fh,
        ):
            cfg = yaml.safe_load(fh)
            cfg["ONT_manifest"] = str(input.manifest)
            cfg["SUNK_len"] = params.sunk_len

            yaml.dump(cfg, out_fh)


rule all:
    input:
        rules.generate_manifest.output,
        rules.generate_config.output,
    default_target:
        True
