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
	params.file_pattern = "/**.bam"
}

if (!params.output_dir) {
	params.output_dir = "gffquant_out"
}
output_dir = "${params.output_dir}/${params.ambig_mode}_${params.mode}"

Channel
	//.fromPath(params.input_dir + "/**.bam")
	//.fromPath(params.input_dir + "/*.fr12RepContigs.sorted.bam")
	.fromPath(params.input_dir + "/" + params.file_pattern)
	.map { file -> 
		def sample = file.name.replaceAll(/\..+/, "")
		return tuple(sample, file)
	}
	.groupTuple()
	.set { samples_ch }

process run_gffquant {
	//conda ""
	publishDir "$output_dir", mode: "link"

	input:
	set sample, file(bamfile) from samples_ch
	//file("${params.db}")

	output:
	stdout result_run_gffquant
	file("${sample}/${sample}.seqname.uniq.txt")
	file("${sample}/${sample}.seqname.dist1.txt")
	file("${sample}/${sample}.feature_counts.txt")
	file("logs/${sample}.o")
	file("logs/${sample}.e")
	
	script:
	"""
	echo $sample $bamfile
	mkdir -p logs
	gffquant ${params.db} ${bamfile} -o ${sample}/${sample} -m ${params.mode} --ambig_mode ${params.ambig_mode} ${params.strand_specific} > logs/${sample}.o 2> logs/${sample}.e
	touch ${sample}/${sample}.seqname.dist1.txt
	"""
}
