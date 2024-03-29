process BLOBTOOLS_ADD {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(blast_out)
    tuple val(meta), path(diamond_proteome_out)
    tuple val(meta), path(diamond_busco_out)
    tuple val(meta), path(assembly_minimap_bam)
    tuple val(meta), path(hic_minimap_bam)
    tuple val(meta), path(lineage1_full_table_tsv)
//    tuple val(meta), path(lineage2_full_table_tsv)
    tuple val(meta), path(blobtools_folder)

    output:
//    tuple val(meta), path("${meta.id}/*.json"), emit: json
    tuple val(meta), path("${meta.id}"), emit: blobtools_folder
    path  "versions.yml"          , emit: versions
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    singularity exec -B /projects ${params.singularity_cache}/blobtoolkit-blobtools_latest.sif blobtools add \
    	--hits ${blast_out} \
	--hits ${diamond_proteome_out} \
        --hits ${diamond_busco_out} \
	--taxrule bestsumorder \
	--taxdump ${params.Blobtoolkit_db}/taxdump \
	--cov ${assembly_minimap_bam} \
	--cov ${hic_minimap_bam} \
	--busco ${lineage1_full_table_tsv} \
	--link taxon.taxid.ENA="https://www.ebi.ac.uk/ena/data/view/Taxon:${params.taxon_taxid}" \
	--link taxon.name.Wikipedia="https://en.wikipedia.org/wiki/${meta.id}" \
    	${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools add: \$(singularity exec -B /projects ${params.singularity_cache}/blobtoolkit-blobtools_latest.sif blobtools --version | sed -e "s/blobtoolkit v//g")
    END_VERSIONS
    """
}
