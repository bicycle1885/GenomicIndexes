    wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
    cat TAIR10_GFF3_genes.gff | sort -k1,1 -k4,4n | bgzip >TAIR10_GFF3_genes.sorted.gff.bgz
    tabix -p gff TAIR10_GFF3_genes.sorted.gff.bgz
    tabix -p gff -C TAIR10_GFF3_genes.sorted.gff.bgz
