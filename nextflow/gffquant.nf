#!/usr/bin/env nextflow

def helpMessage() {
	log.info """
	BLAH!	
	""".stripIndent()
}

if (params.help) {
	helpMessage()
	exit 0
} else {
	params.help = ""
}

if (params.rna) {
	params.strand_specific = "--strand_specific"
} else {
	params.strand_specific = ""
}

if (!params.ambig_mode) {
	params.ambig_mode = "unique_only"
}

if (!params.mode) {
	params.mode = "genome"
}

if (!params.db) {
	if (params.mode == "genome") {
		params.db = "${params.genome_db}"
	} else if (params.mode == "genes") {
		params.db = "${params.genes_db}"
	}
}

if (!params.file_pattern) {
	params.file_pattern = "**.bam"
}

if (!params.publish_mode) {
	params.publish_mode = "link"
}

if (!params.output_dir) {
	params.output_dir = "gffquant_out"
}
output_dir = "${params.output_dir}/${params.ambig_mode}_${params.mode}"

suffix_pattern = params.file_pattern.replaceAll(/\*\*/, "")

Channel
	.fromPath(params.input_dir + "/" + params.file_pattern)
	.map { file -> 
		def sample = file.name.replaceAll(suffix_pattern, "")
		sample = sample.replaceAll(/\.$/, "")
		return tuple(sample, file)
	}
	.groupTuple()
	.set { samples_ch }

process run_gffquant {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	set sample, file(bamfile) from samples_ch

	output:
	stdout result_run_gffquant
	file("${sample}/${sample}.seqname.uniq.txt")
	file("${sample}/${sample}.seqname.dist1.txt")
	file("${sample}/${sample}.feature_counts.txt")
    file("${sample}/${sample}.gene_counts.txt")
	file("logs/${sample}.o")
	file("logs/${sample}.e")
	
	script:
	"""
	echo $sample $bamfile
	mkdir -p logs
	gffquant ${params.db} ${bamfile} -o ${sample}/${sample} -m ${params.mode} --ambig_mode ${params.ambig_mode} ${params.strand_specific} > logs/${sample}.o 2> logs/${sample}.e
	touch ${sample}/${sample}.seqname.dist1.txt ${sample}/${sample}.seqname.uniq.txt ${sample}/${sample}.feature_counts.txt ${sample}/${sample}.gene_counts.txt
	"""
}
