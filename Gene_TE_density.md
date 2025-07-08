````bash
#!/bin/bash

# Define paths
GENOME_FILE="$(pwd)/hap2.chr.fa"
TE_GFF="$(pwd)/argane.asm.chr.full_mask.gff3"
GENE_GFF="$(pwd)/Sideroxylon_spinosum.gff3"
WINDOW_SIZE=500000  # 1 Mb window

# Create a working directory
WORKDIR="Circos_TE_Genes"
mkdir -p $WORKDIR
cd $WORKDIR

# Step 1: Extract chromosome sizes
echo "Extracting chromosome sizes..."
samtools faidx $GENOME_FILE
cut -f1,2 ${GENOME_FILE}.fai > genome.chrom.sizes



# Step 4: Generate genomic windows (0.5 Mb by default)
echo "Generating genomic windows..."
bedtools makewindows -g genome.chrom.sizes -w $WINDOW_SIZE -w 500000 -s 50000 -i winnum > genome_windows.bed

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

karyotype = karyotype.txt  # Path to the karyotype file

<ideogram>
    <spacing>
        default = 0.005r
    </spacing>
    thickness = 60p
    stroke_thickness = 2p
    stroke_color = black
    radius = 0.9r
    show_label     = yes
    label_font     = default
    label_radius   = (dims(ideogram,radius_inner) + dims(ideogram,radius_outer))/2
    label_center   = yes
    label_size     = 20
    label_with_tag = yes
    label_parallel = no
    label_case     = upper


</ideogram>

<plots>
    <plot>
        type = histogram
        file = gene_density.txt # Path to TE density file
        color = green
        thickness = 2p
        r1 = 0.96r
        r0 = 0.76r
        fill_color = green
    </plot>
    <plot>
        type = histogram
        file = te_density.txt # Path to TE density file
        color = red
        thickness = 2p
        r1 = 0.75r
        r0 = 0.55r
        fill_color = red
    </plot>

</plots>



<ticks>
show_ticks       = yes
show_tick_labels = yes

<tick>
spacing        = 50000000
size           = 15p
color          = black
show_label     = yes
label_size     = 24p
label_offset   = 10p
format         = %s
</tick>

</ticks>




















<image>

    <<include /etc/circos/image.conf>>  # Standard image parameters
#    radius           = 1000p    # You can decrease this value (e.g., 500p)
#    image_map_use    = no
#    image_map_overlay = no
#    svg_font_scale   = 0.5      # Decrease font size to match the smaller image size


</image>
<<include /home/slimane/miniconda/pkgs/circos-0.69.9-hdfd78af_0/etc/colors_fonts_patterns.conf>>
<<include /home/slimane/miniconda/pkgs/circos-0.69.9-hdfd78af_0/etc/housekeeping.conf>>

EOL

# Step 8: Run Circos
echo "Running Circos..."
circos -conf circos.conf

echo "Circos plot generated!"



````
