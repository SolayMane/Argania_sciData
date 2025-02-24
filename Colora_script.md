# After preparing the config file, the assembly was proceces usiung this command:
````bash
snakemake --configfile config_argania.yaml --software-deployment-method conda \
--snakefile workflow/Snakefile --cores 56  --rerun-incomplete
````
