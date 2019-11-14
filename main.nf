/*
 * Define the default parameters 
 */
params.outdir     = "/data/results"
params.genome     = "/data/lambda_ref.fasta"
params.reads      = "/data/sample1.fastq"
params.fast5      = "/data/sample1" 

/*
 *  Parse the input parameters
 */

genome_file     = file(params.genome)
reads_ch     = Channel.fromPath(params.reads).map{ it -> [ it.baseName, it  ] }
fast5_ch     = Channel.fromPath(params.fast5).map{ it -> [ it.name, it  ] }
reads_ch.into { nano_reads_ch; filt_reads_ch }


/*
 * Process 1: Nanoplot fastq quality
 */

process 'nanoplot_fastq_quality' { 
  publishDir "$params.outdir/${task.process}_${task.index}/", mode: 'copy', saveAs: { filename -> "report_${name}.html" }

  input:
      tuple val (name), file(reads) from nano_reads_ch

  output:
     tuple val(name), file("NanoPlot-report.html") into nano_out_ch

  script:
  """
  NanoPlot --fastq ${reads} 
  """
}

/*
 * Process 2: Filter fastq 
 */

process 'filtlong_fastq' {
  publishDir 'results_filtlong/', mode: 'copy', saveAs: { filename -> "${name}.filtered.fq" }

  input:
      tuple val(name), file(reads) from filt_reads_ch

  output:
     tuple val(name), file("${name}.filtered.fq") into filtlong_out_ch

  script:
  """
  filtlong --min_mean_q 90 --min_length 3000 ${reads} > sample1.filtered.fq
  """
}

filtlong_out_ch.into { filtlong_denovo_ch; map_filter_ch ; filter_index_ch; call_variants_ch; map_filter_meth_ch}

/*
 * Process 3: Denovo pormoxis
 */

process 'denovo_pomoxis' {
  publishDir 'results_denovo_pomoxis/', mode: 'copy', saveAs: { filename -> "${name}_denovo_pomoxis_final.fa" }

  input:
      tuple val(name), file(filter_reads)  from filtlong_denovo_ch
      

  output:
     tuple val(name), file("denovo/denovo_pomoxis_final.fa") into denovo_pomoxis_out_ch

  script:
  """
  mini_assemble -i ${filter_reads} -o denovo -p denovo_pomoxis -t 8 -c
  """
}


/*
 * Process 4: Minalign reads to reference
 */

process 'minialign_var' {
  publishDir 'results_minialign_var/', mode: 'copy', saveAs: { filename -> "${name}.bam" }

  input:
      file ref_file from genome_file 
      tuple val(name), file(filter_reads) from map_filter_ch


  output:
     tuple val(name), file("${name}.bam"), file("${name}.bam.bai") into minialign_var_ch

  script:
  """
  mini_align -i ${filter_reads} -r ${ref_file} -P -p ${name} -t 8
  samtools index ${name}.bam
  """
}


/*
 * Process 5: Index filtered reads
 */

process 'index_filtered' {
  publishDir 'results_index_filtered/'//, mode: 'copy', saveAs: { filename -> "${name}*" }

  input:
      tuple val(name), file(filter_reads) from filter_index_ch
      tuple val(name), file(fast5_reads) from  fast5_ch   

  output:
     tuple val(name), file("${filter_reads}*") into index_filtered_out_ch

  script:
  """
  nanopolish index -d ${fast5_reads} ${filter_reads}
  """
}

minialign_var_ch.into{ minialign_call_var_ch; minialign_methyl_ch }

/*
 * Process 6: Call variants on polished assembly
 */

//call_variants_ch.view()
//index_filtered_out_ch.view()
call_variants_ch.join(index_filtered_out_ch, by:0)
                .set { fastq_indexed_ch }


fastq_indexed_ch.into{ fastq_indexed_var_ch; fastq_methyl_ch }

process 'call_variants' {
  publishDir 'results_call_variants/', mode: 'copy', saveAs: { filename -> "${name}_variants.vcf" }

  input:
      tuple val(name), file(filter_reads), file(filter_index) from fastq_indexed_var_ch
      file ref_file from genome_file
      tuple val(bam_name), file(filtered_bam), file(filtered_bai) from minialign_call_var_ch
  
  output:
     tuple val(name), file("${name}_variants.vcf") into variants_file_out_ch
  
  script:
  """
  nanopolish variants --consensus -o ${name}_variants.vcf \
    -r ${filter_reads} -b ${filtered_bam} -g ${ref_file} -q dcm,dam --min-candidate-frequency 0.2 -t 8
  """
}

/*
 * Process 7: call methyl
 */

process 'call_methyl' {
  publishDir 'results_call_methyl/', mode: 'copy', saveAs: { filename -> "${name}_methyl.vcf" }

  input:
      tuple val(name), file(filter_reads), file(filter_index) from fastq_methyl_ch
      file ref_file from genome_file
      tuple val(bam_name), file(filtered_bam), file(filtered_bai) from minialign_methyl_ch

  output:
     tuple val(name), file("${name}_*.tsv") into methyl_file_out_ch

  script:
  """
  nanopolish call-methylation -t 8 -r ${filter_reads} -b ${filtered_bam} -g ${ref_file} -q cpg > ${name}_methyl.tsv
  
  python3 /nanopolish/scripts/calculate_methylation_frequency.py -i ${name}_methyl.tsv > ${name}_frequency.tsv 
  """
}

