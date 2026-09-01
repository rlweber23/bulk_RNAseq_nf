#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//samples_csv_ch = Channel.fromPath( params.file_csv )




process trim_adapters_bbuk {
    storeDir params.trimmed_fastq_dir
  
    publishDir { "${params.outdir}/fastqs/trimmed/" }, mode: 'copy'


    input:
    tuple path(fastq_r1), path(fastq_r2)

    output:
    tuple path("${samplename}_R1_trimmed.fastq.gz"),
          path("${samplename}_R2_trimmed.fastq.gz"),
          emit: trimmed_reads

    path "${samplename}.bbduk_stats.txt",
         emit: bbduk_stats


    script:
    samplename = fastq_r1.name.replaceFirst(/_R1_001\.fastq\.gz$/, '')

    """
    ${params.bbduk} \
        in1=${fastq_r1} \
        in2=${fastq_r2} \
        out1=${samplename}_R1_trimmed.fastq.gz \
        out2=${samplename}_R2_trimmed.fastq.gz \
        ref=${params.adapter_fa_file} \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        qtrim=r \
        trimq=${params.trim_q} \
	threads=${task.cpus} \
        stats=${samplename}.bbduk_stats.txt
    """
}




process fastqc_2 {

    publishDir "${params.outdir}/fastqc_post_trim/", mode: 'copy'

    input:
    tuple path(t_fastq_r1), path(t_fastq_r2)

    output:
    path "*_fastqc.*"
 
   script:
    """
    ${params.FASTQC} \
        --threads=${task.cpus} \
        ${t_fastq_r1} \
        ${t_fastq_r2}
    """
}



process multiqc {

    publishDir "${params.outdir}/multiqc_post_trim/", mode: 'copy'

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc \
        . \
        --outdir .
    """
}




process kallisto_index {

    storeDir params.kallisto_index_dir

    input:
    path transcriptome_fasta

    output:
    path "${params.kallisto_index_name}"

    script:
    """
    kallisto index \
        --index=${params.kallisto_index_name} \
        ${transcriptome_fasta} \
        --threads=${task.cpus} 
    """
}


process kallisto_quant {

    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    tuple path(t_fastq_r1), path(t_fastq_r2)

    path kallisto_index
    path gtf
    path chromosome_lengths

    output:
    tuple val(sample_id), path("${sample_id}"),
        emit: kallisto_results

    script:
    sample_id = t_fastq_r1.name.replaceFirst(/_R1_trimmed\.fastq\.gz$/, '')

    genomebam_opts = params.kallisto_genomebam ?
        "--genomebam --gtf ${gtf} --chromosomes ${chromosome_lengths}" :
        ""

    """
    kallisto quant \
        --index=${kallisto_index} \
        --output-dir=${sample_id} \
        --threads=${task.cpus} \
        ${genomebam_opts} \
        ${t_fastq_r1} \
        ${t_fastq_r2}
    """
}


process genomebam_to_bigwig {

    publishDir "${params.outdir}/bigwigs", mode: 'copy'

    input:
    tuple val(sample_id), val(source), path(bam), path(bai)

    output:
    tuple val(sample_id),
          val(source),
          path("${sample_id}.${source}.bw"),
          emit: bigwigs

    script:
    """
    bamCoverage \
        --bam ${bam} \
        --outFileName ${sample_id}.${source}.bw \
        --numberOfProcessors ${task.cpus} \
        --binSize ${params.bigwig_bin_size} \
        --normalizeUsing CPM
    """
}

process make_tx2gene {

    storeDir params.tx2gene_dir

    input:
    path gtf

    output:
    path "${params.tx2gene_name}"

    script:
    """
    zcat ${gtf} | \
    awk -F '\\t' '
        \$3 == "transcript" {
            match(\$9, /gene_id "([^"]+)"/, g)
            match(\$9, /transcript_id "([^"]+)"/, t)

            if (g[1] != "" && t[1] != "")
                print t[1] "\\t" g[1]
        }
    ' > ${params.tx2gene_name}
    """
}




process tximport_counts {

    tag "${counts_from_abundance}"

    publishDir "${params.outdir}/tximport/", mode: 'copy'

    input:
    val counts_from_abundance
    path kallisto_dirs
    path tx2gene

    output:
    tuple val(counts_from_abundance),
          path("gene_counts_${counts_from_abundance}.tsv"),
          emit: gene_counts

    script:

    def output_file = "gene_counts_${counts_from_abundance}.tsv"

    """
    Rscript ${params.tximport_r} \
        ${tx2gene} \
        ${counts_from_abundance} \
        ${output_file}
    """
}


process star_index {

    storeDir params.star_index_store

    input:
    path genome_fasta
    path gtf

    output:
    path "${params.star_index_name}"

    script:
    """
    mkdir -p ${params.star_index_name}

    if [[ "${genome_fasta}" == *.gz ]]; then
        gzip -dc ${genome_fasta} > genome.fa
    else
        cp ${genome_fasta} genome.fa
    fi

    if [[ "${gtf}" == *.gz ]]; then
        gzip -dc ${gtf} > annotation.gtf
    else
        cp ${gtf} annotation.gtf
    fi


    ${params.STAR} \
        --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir ${params.star_index_name} \
        --genomeFastaFiles genome.fa \
        --sjdbGTFfile annotation.gtf \
        --sjdbOverhang ${params.star_sjdb_overhang}
    """

}



process star_align {

    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple path(t_fastq_r1), path(t_fastq_r2)
    path star_index

    output:
    tuple val(sample_id),
          path("${sample_id}.Aligned.sortedByCoord.out.bam"),
          path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"),
          path("${sample_id}.SJ.out.tab"),
          path("${sample_id}.Log.final.out"),
          emit: star_results

    script:
    sample_id = t_fastq_r1.name.replaceFirst(/_R1_trimmed\.fastq\.gz$/, '')

    """
    ${params.STAR} \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --readFilesIn ${t_fastq_r1} ${t_fastq_r2} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}.

    samtools index \
        ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}





workflow {

    /*
     * ============================================================
     * Reference files
     * ============================================================
     */

    transcriptome_fasta = file(params.transcriptome_fasta)
    genome_fasta        = file(params.genome_fasta)
    gtf                 = file(params.gtf)
    chromosome_lengths  = file(params.chromosome_lengths)


    /*
     * ============================================================
     * Reference preparation
     * ============================================================
     */

    kallisto_index_ch = kallisto_index(
        transcriptome_fasta
    )

    tx2gene = make_tx2gene(
        gtf
    )

    if (params.run_star) {

        star_index_ch = star_index(
            genome_fasta,
            gtf
        )
    }


    /*
     * ============================================================
     * Input samples
     * ============================================================
     */

    sample_ch = channel
        .fromPath(params.file_csv)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                file(row.fastq_r1),
                file(row.fastq_r2)
            )
        }


    /*
     * ============================================================
     * Adapter trimming
     * ============================================================
     */

    trimmed_ch = trim_adapters_bbuk(
        sample_ch
    )


    /*
     * ============================================================
     * Optional FastQC / MultiQC
     * ============================================================
     */

    if (params.run_fastqc) {

        fastqc_ch = fastqc_2(
            trimmed_ch.trimmed_reads
        )

        multiqc(
            fastqc_ch.collect()
        )
    }


    /*
     * ============================================================
     * Kallisto quantification
     * ============================================================
     */

    kallisto_ch = kallisto_quant(
        trimmed_ch.trimmed_reads,
        kallisto_index_ch,
        gtf,
        chromosome_lengths
    )


    /*
     * ============================================================
     * Optional STAR alignment
     * ============================================================
     */

    if (params.run_star) {

        star_ch = star_align(
            trimmed_ch.trimmed_reads,
            star_index_ch
        )

        star_bam_ch = star_ch.star_results.map {
            sample_id, bam, bai, sj, log ->

            tuple(
                sample_id,
                "STAR",
                bam,
                bai
            )
        }
    }


    /*
     * ============================================================
     * Optional Kallisto genome BAM
     * ============================================================
     */

    if (params.kallisto_genomebam) {

        kallisto_bam_ch = kallisto_ch.kallisto_results.map {
            sample_id, quant_dir ->

            tuple(
                sample_id,
                "kallisto",
                quant_dir.resolve("pseudoalignments.bam"),
                quant_dir.resolve("pseudoalignments.bam.bai")
            )
        }
    }


    /*
     * ============================================================
     * BigWig generation
     *
     * Combine whichever genome BAM sources were requested,
     * then invoke genomebam_to_bigwig only once.
     * ============================================================
     */

    bigwig_input_ch = channel.empty()

    if (params.run_star) {
        bigwig_input_ch = bigwig_input_ch.mix(star_bam_ch)
    }

    if (params.kallisto_genomebam) {
        bigwig_input_ch = bigwig_input_ch.mix(kallisto_bam_ch)
    }

    if (params.run_star || params.kallisto_genomebam) {

        genomebam_to_bigwig(
            bigwig_input_ch
        )
    }


    /*
     * ============================================================
     * Gene-level quantification
     * ============================================================
     */

    kallisto_dirs_ch = kallisto_ch.kallisto_results.map {
        sample_id, quant_dir ->

        quant_dir
    }

    /*
     * tximport
     */
   
    def counts_from_abundance_options = params.countsFromAbundance ?: ["no"]
   
    // Allow a single string as well as a list
    if (!(counts_from_abundance_options instanceof List)) {
        counts_from_abundance_options = [counts_from_abundance_options]
    }
    
    def valid_counts_from_abundance = [
        "no",
        "scaledTPM",
        "lengthScaledTPM"
    ]
    
    def invalid_counts_from_abundance = counts_from_abundance_options.findAll {
        !(it in valid_counts_from_abundance)
    }
    
    if (invalid_counts_from_abundance) {
        error """
    Invalid countsFromAbundance option(s):
    ${invalid_counts_from_abundance.join(", ")}
    
    Allowed options are:
    ${valid_counts_from_abundance.join(", ")}
    """
    }
    
    counts_from_abundance_ch = Channel.fromList(
        counts_from_abundance_options
    )
    
    tximport_counts(
        counts_from_abundance_ch,
        kallisto_dirs_ch.collect(),
        tx2gene
    )
   

}

