process BED_PROCESSING {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_sorted.bed"), emit: sorted_bed

    script:
    """
    #awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"/1""\t"\$5"\t"\$6}' $bed_F > F.bed
    #awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"/2""\t"\$5"\t"\$6}' $bed_R > R.bed
    #cat F.bed R.bed > concat.bed
    #for BAM in `ls $projectDir/results/hicpro/mapping/*.bam`;  
			#do bedtools bamtobed -i $BAM > merged_bed

    bedtools bamtobed -i $bam > merged_bed
    sort --parallel=8 --buffer-size=80%  --temporary-directory=$projectDir --output=bed_sorted.bed \$merged_bed
    """
}
