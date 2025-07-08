# Jupiter plot
````bash
JupiterPlot-1.1/jupiter name=SideroVsParadox t=32 fa=vitelariaParadoxa_chr.fasta ref=hap2.chr.fa ng=100 labels=both
````


# After preparing the config file, the assembly was proceces usiung this command:
````bash
snakemake --configfile config_argania.yaml --software-deployment-method conda \
--snakefile workflow/Snakefile --cores 56  --rerun-incomplete
````
