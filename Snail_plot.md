# to generate a snail plot from the assembled genome we used blobtools kit
````bash
blobtools create --fasta yahs_primary.filtered.fasta --busco filtered_busco.tsv BlobDirFiltred
blobtools view --plot --view snail BlobDirFiltred/

````
