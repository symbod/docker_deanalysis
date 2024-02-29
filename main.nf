#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = '' // Path to meta.txt
params.count_file = '' // Path to count.txt

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05

// scripts
deanalysis_script = Channel.fromPath("${projectDir}/DEAnalysis.R")
summarizer_script = Channel.fromPath("${projectDir}/Summarizer.R")


process deanalysis {
    container 'kadam0/deanalysis:0.0.1'
    publishDir params.output, mode: "copy"

    input:
    path script_file from deanalysis_script
    path meta_file from Channel.fromPath(params.meta_file)
    path count_file from Channel.fromPath(params.count_file)

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --count_file ${count_file} --out_dir ${params.output} --logFC ${params.logFC} --logFC_up ${params.logFC_up} --logFC_down ${params.logFC_down} --p_adj ${params.p_adj} --alpha ${params.alpha}
    """
}

process summarize {
    container 'kadam0/deanalysis:0.0.1'
    publishDir params.output, mode: "copy"

    input:
    path script_file from summarize_script
    val input_dir from params.output
    val logFC from params.logFC
    val logFC_up from params.logFC_up
    val logFC_down from params.logFC_down
    val p_adj from params.p_adj
    val alpha from params.alpha

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --out_dir ${input_dir} --logFC ${logFC} --logFC_up ${logFC_up} --logFC_down ${logFC_down} --p_adj ${p_adj} --alpha ${alpha}    """
}



workflow {
  deanalysis()
  summarize()
}
