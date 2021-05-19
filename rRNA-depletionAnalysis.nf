
def helpMessage() {
    log.info"""
    ==------------------------------------==
    Personal use RNA pipeline
    ==------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams  "*.bam" --gff <.gff.gz> -rRNAgenes A,B,C -profile conda

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
params.paired = false
params.featureCountsS = 'reverse'
params.g = 'locus_tag'
params.t = 'transcript'
params.extraAttributes = 'gene_name'
params.pubDir = "Results"
pubDir = file(params.pubDir)
gff_file = file(params.gff)
alignment_files = file(params.bams)

isPaired = params.paired
fcStrandness = params.featureCountsS
featureCounts_g = params.g
featureCounts_t = params.t
featureCounts_extra = params.extraAttributes
rRNA_list = params.rRNAgenes


process unZipGFF{
   input:
      file zipped_gff from gff_file
   output:
      file "${zipped_gff.baseName}" into gff3
   script:
      """
      gunzip -c $zipped_gff > ${zipped_gff.baseName}
      """
}
/*
Convert GFF to GTF
*/
process convertGFFtoGTF {
    tag "$gff"
    input:
    file gff from gff3

    output:
       file "${gff.baseName}.gtf" into gtf
    script:
       """
       gffread $gff -F -T -o ${gff.baseName}.gtf
       """
}
// 4) counts features for depletion analysis
process featureCounts {
  input:
    val g from featureCounts_g
    val t from featureCounts_t
    val extraAttributes from featureCounts_extra
    file annotation from gtf
    file(bams) from alignment_files.collect()
    publishDir "$pubDir/rRNADepletion/", mode: 'copy'
  output:
    file("counts.txt") into feature_counts
    file("counts.txt.summary") into featureCounts_logs
  script:
    s = "0"
    p = ""
    if(fcStrandness == "forward"){
      s = "1"
    } else if(fcStrandness == "reverse"){
      s = "2"
    }
    if(isPaired){
      p ="-p"
    }
    """
    featureCounts -a $annotation -s $s -O -M --fraction -t $t -g $g --extraAttributes $extraAttributes --minOverlap 10 --fracOverlap 0.9 $p $bams
    """
}

process depletionCalculation{
  input:
    file featureCounts from feature_counts
    val g from featureCounts_g
    val list from rRNA_list
    publishDir "$pubDir/rRNADepletion", mode: 'copy'
  output:
    file("rRNA_remaining.csv") into depletion_calc
  script:
  """
    compute_ratio.py -c $featureCounts -r $list
  """
}
