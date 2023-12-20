#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = '' // Path to meta.txt
params.count_file = '' // Path to count.txt

// scripts
deanalysis_script = Channel.fromPath("${projectDir}/DEAnalysis.R")


process deanalysis {
    container 'kadam0/deanalysis:0.0.1'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path meta_file 
    path count_file 

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --count_file ${count_file} --out_dir ./
    """
}



workflow {
  meta_file = Channel.fromPath(params.meta_file)
  count_file = Channel.fromPath(params.count_file)

  deanalysis(deanalysis_script, meta_file, count_file)
}
