#!/bin/
def helpMessage() {
    log.info"""
    ==------------------------------------==
    Personal use RNA pipeline
    ==------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run rRNA-depletionAnalysis.nf --bams  "*.bam" --gff <.gff.gz> -rRNAgenes A,B,C -profile conda

    Mandatory arguments:
      --bams                  .bam input files containg the reads
      --rRNAgenes             Input String of rRNA genes separated by ","
      --gff                   Input file in gff (.gz) format, annotations

    Optional arguments:
      --paired                Paired end reads
      --featureCountsS        FeatureCounts attribute strandedness: unstranded, reverse or forward strandness (default: reverse)
      --g                     FeatureCounts attribute to group features (default: locus_tag)
      --t                     FeatureCounts attribute that should be counted (default: transcript)
      --extraAttributes       FeatureCounts extra attriutes divided by comma (default: gene_name)
      --pubDir                The directory where the results will be stored [def: Results]
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}
params.paired = false
params.featureCountsS = 'reverse'
params.g = 'locus_tag'
params.t = 'transcript'
params.extraAttributes = 'gene_name'
params.pubDir = "Results"
pubDir = file(params.pubDir)
gff_file = Channel.fromPath(params.gff, checkIfExists: true)
alignment_files = Channel.fromPath(params.bams, checkIfExists:true)

isPaired = params.paired
fcStrandness = params.featureCountsS
featureCounts_g = params.g
featureCounts_t = params.t
featureCounts_extra = params.extraAttributes
rRNA_path = Channel.fromPath(params.rRNAgenes, checkIfExists:true)
 
process unZipGFF{
    debug true

    input:
    path zipped_gff

    output:
    path "${zipped_gff.baseName}"

    script:
    """
    gunzip -c $zipped_gff > ${zipped_gff.baseName}
    """
}

// Convert GFF TO GTF
process convertGFFtoGTF{
    conda 'environment.yml'
    tag "$gff"

    input:
    path gff

    output:
    path "${gff.baseName}.gtf"

    script:
    """
    gffread $gff -F -T -o ${gff.baseName}.gtf
    """

}

// counts features for depletion analysis
process featureCounts{
    publishDir "./rRNADepletion", mode: 'copy'
    conda 'environment.yml'

    input:
    val(g)
    val(t)
    val(extraAttributes)
    path(annotation)
    path(bams)

    output:
    tuple path("counts.txt"), path("counts.txt.summary")

    script:
    s = "0"
    p = ""

    if(fcStrandness == "forward"){
        s  = "1"
    } else if (fcStrandness == "reverse"){
        s = "2"
    }

    if (isPaired){
        p = "-p"
    }

    """
    featureCounts -a $annotation -s $s -O -M --fraction -t $t -g $g --extraAttributes $extraAttributes -o counts.txt --minOverlap 10 --fracOverlap 0.9 $p $bams
    """
}

process depletionCalculation{
    debug true
    conda 'environment.yml'
    publishDir "./rRNADepletion", mode: 'copy'

    input:
    path(featureCounts)
    val g
    path(list) 

    output:
    path("rRNA_remaining.csv")

    script:
    """
    compute_ratio.py -c $featureCounts -r $list
    """
}

process tpmCalculation{
    conda 'environment.yml'
    publishDir "./rRNADepletion", mode: 'copy'

    input:
    path(bamFile)
    path(gffFile)
    path(gtfFile)

    output:
    path '*.results' 

    script:
    """
    rsem-prepare-reference --gtf $gtfFile 
    rsem-calculate-expression -p 8 --bam --estimate-rspd --append-names --output-genome-bam $bamFile $gffFile ./rRNADepletion
    """
}


workflow {  
    ch_gff = unZipGFF(gff_file)
    ch_gtf = convertGFFtoGTF(ch_gff)
    
    ch_feature_counts = featureCounts(featureCounts_g, featureCounts_t, featureCounts_extra, ch_gtf, alignment_files.collect())
    ch_feature_counts
        | flatten
        | branch { it ->
            txt: it.fileName.endsWith('counts.txt')
            summary: it.fileName.endsWith('counts.txt.summary')
        }
        | set { result }

    ch_depletion_calculation = depletionCalculation(result.txt, featureCounts_g, rRNA_path)
    //isoform_results = tpmCalculation(alignment_files, ch_gff, ch_gtf)
}