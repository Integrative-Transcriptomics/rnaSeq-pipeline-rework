
def helpMessage() {
    log.info"""
    ==------------------------------------==
    Personal use RNA pipeline
    ==------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --reads  "*.fastq.gz" --reference <.fna.gz> --gff <.gff.gz> -profile conda

    Mandatory arguments:
      --reads                 FastQ input file containg the reads
      --reference             Input File in fasta format (gz), nucleotide sequences. Each entry (">") is considered as one gene.
      --gff                   Input file in gff format, annotations

    Optional arguments:
      --paired                Paired end reads
      --noQuali               Disable quality control
      --noCounts              Disable feature counts
      --featureCountsS        unstranded, reverse or forward strandness (default: reverse)
      --g                     FeatureCounts attribute to group features (default: locus_tag)
      --t                     FeatureCounts attribute that should be counted (default: transcript)
      --extraAttributes       FeatureCounts extra attriutes divided by comma (default: gene_name)
      --M                     Count multi-mapping reads
      --O                     Count reads overlapping features
      --fraction              Fractional counts for multi-mapping/overlapping features (must be used together with -M or -O)
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
params.noQuali = false
params.noCounts = false
params.featureCountsS = 'reverse'
params.g = 'locus_tag'
params.t = 'transcript'
params.extraAttributes = 'gene_name'
params.M = false
params.O = false
params.fraction = false
if(params.fraction && !(params.O || params.M)){
  helpMessage()
  exit 0
}
params.pubDir = "Results"
pubDir = file(params.pubDir)
genome_file = file(params.reference)
gff_file = file(params.gff)

isPaired = params.paired
noQualiControl = params.noQuali
noFC = params.noCounts
fcStrandness = params.featureCountsS
featureCounts_g = params.g
featureCounts_t = params.t
featureCounts_extra = params.extraAttributes
multiMapping = params.M
overlapping = params.O
fraction = params.fraction

// The basic input of the Pipeline. based on the Read name, a <base>_reference.info file is locatad, containing the information of the reference to be used for this sample
Channel.fromFilePairs(params.reads, size: 1)
        .ifEmpty { exit 1, "Readfiles not specified" }
        .into { reads_fastQC; reads_trimgalore }

/*
Unzip genome file
*/
process unZipGenome {
  input:
    file zipped_genome from genome_file
  output:
    file "${zipped_genome.baseName}" into genome_fasta
  script:
    """
      gunzip -c $zipped_genome > ${zipped_genome.baseName}
    """
}

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
/*
* FastQC for raw reads
*/
process fastqc {
  tag "$id"
  label 'process_medium'
  publishDir "$pubDir/fastqc", mode: 'copy',
     saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }
  when:
  noQualiControl == false
  input:
    tuple id, file(reads) from reads_fastQC

  output:
    file "*_fastqc.{zip,html}" into fastqc_results

  script:
    """
      fastqc --threads ${task.cpus} --quiet $reads
    """
}


// 1) Trimm the reads
process trimming {
  tag "$id"
  input:
    tuple id, file( read_file ) from reads_trimgalore
    publishDir "$pubDir/Trimgalore/", mode: 'copy',
      saveAs: {filename ->
        if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
        else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
        else null
}

  output:
    set id, file("*_trimmed.fq") into ch_trimmed
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_results
    file "*trimming_report.txt" into trimgalore_results

  script:
  if(!isPaired){
    """
      trim_galore --dont_gzip --fastqc $read_file
    """
  }else{
    """
      trim_galore --dont_gzip --paired --retain_unpaired --fastqc ${read_file[0]} ${read_file[1]}
    """
  }
}

// 2) Create hisat_index
process hisat2_index {
  input:
    file reference_genome from genome_fasta
  output:
    file "${reference_genome.baseName}.*.ht2*" into hisat2_indeces
  script:
    """
      hisat2-build -p ${task.cpus} $reference_genome ${reference_genome.baseName}.hisat2_index
    """
}

// 3) Create hisat_mapping
process hisat2_mapping {
    input:
      file indeces from hisat2_indeces.collect()
      tuple id, file(reads) from ch_trimmed
    output:
      tuple id, file("*.bam")into alignment_files
      file("*.bam")into alignment_files_fc
      tuple id, file("*_summary.txt") into alignment_logs
    script:
      index_name = indeces[0].toString() - ~/.\d.ht2l?/
      if(!isPaired){
      """
      hisat2 -x $index_name \
                   -U $reads \
                   --no-spliced-alignment \
                   -p ${task.cpus} \
                   --met-stderr \
                   --new-summary \
                   --summary-file ${id}_summary.txt --rg-id ${id} --rg SM:${id} \
                   | samtools view -bS -F 4 -F 8 -F 256 - > ${id}.bam
      """
    } else{
      """
      hisat2 -x $index_name \
                   --1 ${reads[0]} \
                   --2 ${reads[1]} \\
                   -p ${task.cpus} \
                   --met-stderr \
                   --new-summary \
                   --summary-file ${id}_summary.txt --rg-id ${id} --rg SM:${id} \
                   | samtools view -bS -F 4 -F 8 -F 256 - > ${id}.bam
      """
    }
}

// 4) counts features
process featureCounts {
  input:
    val g from featureCounts_g
    val t from featureCounts_t
    val extraAttributes from featureCounts_extra
    file annotation from gtf
    file(bams) from alignment_files_fc.collect()
    publishDir "$pubDir/FeatureCounts/", mode: 'copy'
    when:
    noFC == false
  output:
    file("counts.txt") into feature_counts
    file("counts.txt.summary") into featureCounts_logs
  script:
    M = ""
    O = ""
    frac = ""
    s = "0"
    if(fcStrandness == "forward"){
      s = "1"
    } else if(fcStrandness == "reverse"){
      s = "2"
    }
    if(multiMapping){
      M = "-M"
    }
    if(overlapping){
      O = "-O"
    }
    if(fraction){
      frac = "--fraction"
    }
    """
    featureCounts -a $annotation -s $s -t $t -g $g --extraAttributes $extraAttributes -o counts.txt $M $O $frac $bams
    """
}

// samtools sort index
process samtools {
  input:
    tag "$id"
    cpus 4
    publishDir "$pubDir/Bams/", mode: 'copy'

  input:
    tuple id, file(bam) from alignment_files

  output:
    tuple id, file ("${id}.sorted.bam*") into sorted_alignment_files_bamqc
    tuple id, file ("${id}.sorted.bam*") into sorted_alignment_files_rnaseq


  script:
  """
    samtools sort -@ ${task.cpus} -o ${id}.sorted.bam $bam
    samtools index ${id}.sorted.bam
  """
}

// Qualimap auf Bam
process qualimap_bamqc {
  tag "$id"
  publishDir "$pubDir/QualiMapsBamQC/", mode: 'copy'
  when:
  noQualiControl == false
  input:
    tuple id, file(bam) from sorted_alignment_files_bamqc

  output:
    tuple id, file("${id}") into qualimap_bamqc_results

  script:
  sorted_bam = "${id}.sorted.bam"
  """
    qualimap bamqc -nt ${task.cpus} -bam $sorted_bam -outdir ${id}
  """
}

process qualimap_rnaseq {
  tag "$id"
  publishDir "$pubDir/QualiMapsRNAseq/", mode: 'copy'
  when:
  noQualiControl == false
  input:
    tuple id, file(bam) from sorted_alignment_files_rnaseq
    file annotation from gtf
  output:
    tuple id, file("${id}") into qualimap_rnaseq_results
  script:
  sorted_bam = "${id}.sorted.bam"
    """
      qualimap rnaseq -bam $sorted_bam -gtf $annotation -outdir ${id}
    """
}

// MultiQC
process multiqc {
    publishDir "${pubDir}/MultiQC", mode: 'copy'

    when:
    noQualiControl == false
    input:
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
    file ('trimgaloreFastQC/*') from trimgalore_fastqc_results.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file ('qualimapBAMQC/*') from qualimap_bamqc_results.collect().ifEmpty([])
    file ('qualimapRNAseq/*') from qualimap_rnaseq_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into multiqc_report

    script:
    """
    multiqc .\\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m sortmerna -m fastqc -m qualimap -m salmon
    """
}
