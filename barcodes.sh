// link config file with '-c config_file'

bc_script='/ei/workarea/group-scg/wojtowic/SS2_26032020/barcodes/extract_barcode_from_bbduk_nnn.R' //uses tidyverse
aggregating_script='/ei/workarea/group-scg/wojtowic/SS2_26032020/barcodes/aggregating_barcode_lists.R' //uses tidyverse
image_scater='/ei/workarea/group-ga/uzun/scqc/Containers/R-3.5.2_scater.img'


Channel
	.fromFilePairs("${params.reads}${params.pattern}.fastq.gz")
	.set{read_pairs}
	
process bbduk_filter {
    publishDir "$params.barcodes_dir"
    beforeScript 'source bbmap-37.24'
    tag "$sampleID"

    input:
    set val(sampleId), file(reads) from read_pairs
	
    output:
    set val(sampleId), file("match_${sampleId}_bb.fq") into matched_seqs

    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} out=nonmatch_${sampleId}_bb.fq outm=match_${sampleId}_bb.fq literal=${params.literal}  k=${params.literal_length} hdist=1 mm=f rcomp=t threads=1

    """
}

process bbduk_mask {

	publishDir "$params.barcodes_dir"
	beforeScript 'source bbmap-37.24'
	
	input:
	set sampleId, file(matched_seq) from matched_seqs

	output:
	set file("mask_proc_${sampleId}.match.fq"), file("mask_proc_${sampleId}.nonm.fq") into masked

	"""
	bbduk.sh in=${matched_seq} out=mask_proc_${sampleId}.nonm.fq outm=mask_proc_${sampleId}.match.fq literal=${params.literal} kmask=N  k=${params.literal_length} hdist=1 mm=f rcomp=t threads=1
	"""
}

process r_prep {

input:
file(masked_bbs) from masked.filter{ it[1].size()>0 }

output:

file("r_prepped_${masked_bbs[1]}.txt") into r_prepped

"""
cat ${masked_bbs[1]} | grep -E '^[TCGAN]+' > r_prepped_${masked_bbs[1]}.txt
"""
}

process bc_extract {

publishDir "$params.final_path"

	input:
	file r_prep_file from r_prepped

	output:

	file("final_bc_${r_prep_file}") into last_channel


"""
singularity exec ${image_scater} Rscript ${bc_script} ${r_prep_file} final_bc_${r_prep_file} ${params.literal_length}
"""

}

process bc_aggregate {

publishDir "$params.aggregated_file_dir"

	input:
	file flag_check from last_channel.collect()

"""
singularity exec ${image_scater} Rscript ${aggregating_script} ${params.final_path} ${params.agg_filename}
"""
}
