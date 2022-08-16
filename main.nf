#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def classify_sample(sample, files) {

    def meta = [:]
    meta.is_paired = (files instanceof Collection && files.size() == 2)
    meta.id = sample

    return [meta, files]

    if (meta.is_paired) {
        return [meta, files]
    }

    return [meta, [files]]

}


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
	} else if (params.mode == "domain") {
		params.db = "${params.domain_db}"
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


process stream_minimap2_gffquant {
	publishDir "${output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(fastq)
	path(index)
	path(db)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_params = "-o ${sample.id}/${sample.id} -m ${params.mode} --ambig_mode ${params.ambig_mode} ${params.strand_specific}"
	def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"
	def mm_options = "--sam-hit-only -t ${task.cpus} -x sr --secondary=yes -a"
	"""
	mkdir -p logs/
	echo 'Copying database...'
	cp -v ${db} gq_db.sqlite3
	minimap2 ${mm_options} --split-prefix ${sample.id}_split ${index} ${reads} | gffquant ${gq_params} gq_db.sqlite3 - > logs/${sample}.o 2> logs/${sample}.e
	rm -v gq_db.sqlite3
	"""
	// minimap2 --sam-hit-only -t <threads> -x sr --secondary=yes -a [-o <out.sam>] --split-prefix <prefix> <mmi> <reads>
}


process run_gffquant {
	publishDir "${output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(bam)
	path(db)

	output:
	tuple val(sample), path("${sample}/*.txt.gz"), emit: results

	script:
	def gq_params = "-o ${sample}/${sample} -m ${params.mode} --ambig_mode ${params.ambig_mode} ${params.strand_specific}"
	gq_params += (params.unmarked_orphans) ? " --unmarked_orphans" : ""
	if (params.do_name_sort) {
		"""
		mkdir -p logs/
		echo 'Copying database...'
		cp -v ${db} gq_db.sqlite3
		samtools collate -O ${bam} -@ ${task.cpus} | gffquant ${gq_params} gq_db.sqlite3 - > logs/${sample}.o 2> logs/${sample}.e
		rm -v gq_db.sqlite3
		"""
	} else {
		"""
		mkdir -p logs/
		echo 'Copying database...'
		cp -v ${db} gq_db.sqlite3
		gffquant ${gq_params} gq_db.sqlite3 ${bam} > logs/${sample}.o 2> logs/${sample}.e
		rm -v gq_db.sqlite3
		"""
	}
}

process collate_feature_counts {
	publishDir "${output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(count_tables)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/
	collate_counts . -o collated/collated -c uniq_scaled
	collate_counts . -o collated/collated -c combined_scaled
	"""
}


workflow {

	feature_count_ch = Channel.empty()

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + params.file_pattern)
		.map { file ->
			def sample = file.name.replaceAll(suffix_pattern, "")
			sample = sample.replaceAll(/\.$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort:true)

	run_gffquant(bam_ch, params.db)
	feature_count_ch = feature_count_ch
		.concat(run_gffquant.out.results)

	fastq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{fastq.gz,fq.gz}")
		.map { file ->
			def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
			sample = sample.replaceAll(/_R?[12]$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)
		.map { classify_sample(it[0], it[1]) }

	if (params.minimap2_index) {
		stream_minimap2_gffquant(fastq_ch, params.minimap2_index, params.db)
		feature_count_ch = feature_count_ch
			.concat(stream_minimap2_gffquant.out.results)
	}

	feature_count_ch = feature_count_ch
		.map { sample, files -> return files }
		.flatten()
		.filter { !it.name.endsWith("gene_counts.txt") }
		.filter { !it.name.endsWith("seqname.uniq.txt") }
		.filter { !it.name.endsWith("seqname.dist1.txt") }
		.map { file -> 
			def category = file.name.replaceAll(/\.txt$/, "")
				.replaceAll(/.+\./, "")
			return tuple(category, file)
		}
		.groupTuple(sort:true)

	//feature_count_ch.view()

	if (!params.no_collate) {

		collate_feature_counts(feature_count_ch)

	}

}
