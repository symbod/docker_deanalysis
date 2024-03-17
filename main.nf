#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = '' // Path to meta.txt
params.count_file = '' // Path to count.txt

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05

// data
meta_file = file(params.meta_file)
count_file = file(params.count_file)

// scripts
deanalysis_script = file("${projectDir}/DEAnalysis.R")
summarizer_script = file("${projectDir}/Summarizer.R")

process deanalysis {
    container 'kadam0/deanalysis:0.0.2'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path meta_file 
    path count_file 

    output:                                
    path("de_out", type:"dir")

    script:
    """
    Rscript $script_file --meta_file $meta_file --count_file $count_file --out_dir ./de_out --logFC $params.logFC --logFC_up $params.logFC_up --logFC_down $params.logFC_down --p_adj $params.p_adj --alpha $params.alpha
    """
}

process summarize {
    container 'kadam0/deanalysis:0.0.2'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path de_output

    output:                                
    path "summary"

    script:
    """
    Rscript $script_file --in_dir $de_output --out_dir ./summary --logFC $params.logFC --logFC_up $params.logFC_up --logFC_down $params.logFC_down --p_adj $params.p_adj --alpha $params.alpha
    """
}

workflow {
  deanalysis(deanalysis_script, meta_file, count_file)
  summarize(summarizer_script, deanalysis.out)
}
