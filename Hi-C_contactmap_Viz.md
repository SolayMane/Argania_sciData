# How to generate hic.out
<code>juicer pre results/yahs_primary/asm_yahs.bin results/yahs_primary/asm_yahs_scaffolds_final.agp \
results/bwa_index_primary/asm.fa.fai | \
sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part" </code>

<code> java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt.part out.hic.part scaffolds_final.chrom.sizes && mv out.hic.part out.hic </code>
# Finally, the output file out.hic could be used for visualisation with Juicebox.
