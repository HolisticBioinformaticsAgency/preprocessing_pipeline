process FASTQC {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

  // publish per subject
  publishDir { "${params.outdir}/${subject}/fastqc" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(reads)   // reads = [R1,R2] (paired)

  output:
  path "*_fastqc.html", emit: html
  path "*_fastqc.zip",  emit: zip

  script:
  """
  fastqc -q -o ./ ${reads.join(' ')}
  """
}
