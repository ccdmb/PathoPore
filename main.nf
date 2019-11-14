/*
 * Define the default parameters 
 */

params.genome     = "/data/lambda_ref.fasta"
params.reads      = "/data/pass.fastq"
params.fast5      = "/data/fast5_reads" 

/*
 *  Parse the input parameters
 */

genome_file     = file(params.genome)
reads_ch     = Channel.fromPath(params.reads).map{ it -> [ it, it.name  ] }
fast5_ch     = Channel.fromPath(params.fast5).map{ it -> [ it, it.name  ] }
reads_ch.into { nano_reads_ch; filt_reads_ch }


/*
 * Process 1: Nanoplot fastq quality
 */

process 'nanoplot_fastq_quality' { 
  publishDir 'results_nanoplot/', mode: 'copy', saveAs: { filename -> "report_${name}.html" }

  input:
      tuple file(reads), val(name) from nano_reads_ch

  output:
     tuple file("NanoPlot-report.html"), val(name) into nano_out_ch

  script:
  """
  NanoPlot --fastq ${reads} 
  """
}

/*
 * Process 2: Filter fastq 
 */

process 'filtlong_fastq' {
  publishDir 'results_filtlong/', mode: 'copy', saveAs: { filename -> "${name}_filter.fq" }

  input:
      tuple file(reads), val(name) from filt_reads_ch

  output:
     tuple file("filter.fq"), val(name) into filtlong_out_ch

  script:
  """
  filtlong --min_mean_q 90 --min_length 3000 ${reads} > filter.fq
  """
}

filtlong_out_ch.into { filtlong_denovo_ch; map_filter_ch ; filter_index_ch}

/*
 * Process 3: Denovo pormoxis
 */

process 'denovo_pomoxis' {
  publishDir 'results_denovo_pomoxis/', mode: 'copy', saveAs: { filename -> "${name}_denovo_pomoxis_final.fa" }

  input:
      tuple file(filter_reads), val(name) from filtlong_denovo_ch
      

  output:
     tuple file("denovo/denovo_pomoxis_final.fa"), val(name) into denovo_pomoxis_out_ch

  script:
  """
  mini_assemble -i ${filter_reads} -o denovo -p denovo_pomoxis -t 8 -c
  """
}


/*
 * Process 4: Minalign reads to reference
 */

process 'minialign_var' {
  publishDir 'results_minialign_var/', mode: 'copy', saveAs: { filename -> "${name}_minialign_var.bam" }

  input:
      file ref_file from genome_file 
      tuple file(filter_reads), val(name) from map_filter_ch


  output:
     tuple file("minialign_var.bam"), val(name) into minialign_var_ch

  script:
  """
  mini_align -i ${filter_reads} -r ${ref_file} -P -p minialign_var -t 8
  """
}


/*
 * Process 5: Index filtered reads
 */

process 'index_filtered' {
  publishDir 'results_index_filtered/'

  input:
      tuple file(filter_reads), val(name) from filter_index_ch
      tuple file(fast5_reads), val(name) from  fast5_ch   

  output:
     tuple file(fast5_reads), val(name) into index_filtered_out_ch

  script:
  """
  nanopolish index -d ${fast5_reads} ${filter_reads}
  """
}

