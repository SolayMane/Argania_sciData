````bash
#!/bin/bash

# Define paths
GENOME_FILE="/ArganiaGenomics/Assembly_Colora/results/assemblies/yahs_primary.fa"
TE_GFF="/ArganiaGenomics/Assembly_Colora/TE/RepeatMasker_out/yahs_primary.fa.out.gff"
GENE_GFF="/ArganiaGenomics/Plots/te_gene/yahs_primary_genes_features.gff3"
WINDOW_SIZE=1000000  # 1 Mb window

# Create a working directory
WORKDIR="/ArganiaGenomics/Plots/te_gene"
mkdir -p $WORKDIR
cd $WORKDIR

# Step 1: Extract chromosome sizes
echo "Extracting chromosome sizes..."
samtools faidx $GENOME_FILE
cut -f1,2 ${GENOME_FILE}.fai > genome.chrom.sizes



# Step 4: Generate genomic windows (1 Mb by default)
echo "Generating genomic windows..."
bedtools makewindows -g genome.chrom.sizes -w $WINDOW_SIZE -s 100000 -i winnum > genome_windows.bed

# Step 5: Calculate TE density in each window
echo "Calculating TE density..."
bedtools coverage -a genome_windows.bed -b $TE_GFF > te_density.bed

# Step 6: Calculate Gene density in each window
echo "Calculating Gene density..."
bedtools coverage -a genome_windows.bed -b $GENE_GFF > gene_density.bed

# Step 7: Prepare TE density for Circos
echo "Preparing TE density for Circos..."
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $8}' te_density.bed > te_density.txt

# Step 8: Prepare Gene density for Circos
echo "Preparing Gene density for Circos..."
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $8}' gene_density.bed > gene_density.txt

# Step 9: Create karyotype file for Circos
echo "Creating karyotype file..."
awk 'BEGIN {OFS="\t"} {print "chr", "-", $1, $1, 0, $2, "chr"$1}' genome.chrom.sizes > karyotype.txt


cat <<EOL > circos.conf
# circos.conf
karyotype = karyotype_clean.txt   # Path to the karyotype file
chromosomes_units  = 1000000
chromosomes_display_default = yes  # Do not show labels for all chromosomes


<ideogram>
    <spacing>
        default = 0.003r
    </spacing>
    thickness = 60p
    stroke_thickness = 2p
    stroke_color = black
    radius = 0.9r

    show_label     = yes
    label_font     = bold
    label_radius   = (dims(ideogram,radius_inner) + dims(ideogram,radius_outer))/1.8
    label_center   = yes
    label_size     = 35
    label_with_tag = yes
    label_parallel = yes
    label_case     = lower
    #label_format     = eval( replace(var(label),"scaffold_","Chr") )
#    label_format     = eval( var(chr) =~ /scaffold_[1-11]$/ ? var(label) : "")
    label_format = eval(var(chr) =~ /scaffold_(?:[1-9]|10)$/ ? replace(var(label),"scaffold_","Chr") : "")

</ideogram>

show_ticks = yes
show_tick_labels = yes

<ticks>
    radius           = 1r
    color            = black
    thickness        = 2p
    multiplier       = 1e-6  # Convert base pairs to Mb
    format           = %d Mb
    <tick>
        spacing     = 10u  # Every 10 Mb
        size        = 15p
        thickness   = 2p
        color       = black
        label_offset = 2p
        show_label  = yes
        label_size  = 15p
    </tick>
</ticks>
<plots>
    <plot>
        type = histogram
        file = gene_density_clean.txt # Path to TE density file
        color = green
        thickness = 2p
        r1 = 0.97r
        r0 = 0.9r
        fill_color = green
    </plot>
    <plot>
        type = histogram
        file = te_density_clean.txt # Path to TE density file
        color = red
        thickness = 2p
        r1 = 0.89r
        r0 = 0.82r
        fill_color = red
    </plot>
    <plot>
        type = line
        file = gc_content.txt
        color = dark
        thickness = 2p
        r1 = 0.7r
        r0 = 0.6r
        fill_color = dark
#        max = 0.60
#       min = 0.20
    </plot>
</plots>

<image>
    <<include /etc/circos/image.conf>>  # Standard image parameters
</image>

<<include /etc/circos/colors_fonts_patterns.conf>>
<<include /etc/circos/housekeeping.conf>>



````
