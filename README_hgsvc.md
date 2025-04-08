To generate config and split reads/assembly.
```bash
snakemake -np -s generate_config.smk --config run_name="" sunk_len=30
```

Run GAVISUNK.
```bash
snakemake -np --workflow-profile workflow/profiles/lpc/ --configfile config_hgsvc.yaml -s workflow/Snakefile
```