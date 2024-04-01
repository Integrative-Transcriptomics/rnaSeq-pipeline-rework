#!/bin/
def helpMessage() {
    log.info"""
    ==------------------------------------==
    Personal use RNA pipeline
    ==------------------------------------==
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --reads  "*.fastq.gz" --reference <.fna.gz> --gff <.gff.gz> -profile conda

    Mandatory arguments:
      --reads                 FastQ (.gz) input file containg the reads
      --reference             Input File in fasta format (.gz), nucleotide sequences.
      --gff                   Input file in gff (.gz) format, annotations
      --rRNAgenes             Input file in txt format, header: gene,type; body: rRNA genes (example file: rRNA_genes_type.txt)

    Optional arguments:
      --paired                Paired end reads
      --noQuali               Disable quality control
      --noQualiRNA            Disable quality control RNASeq
      --noCounts              Disable feature counts
      --noDepletionQC           Disable depletion Analysis
      --featureCountsS        FeatureCounts attribute strandedness: unstranded, reverse or forward strandness (default: reverse)
      --g                     FeatureCounts attribute to group features (default: locus_tag)
      --t                     FeatureCounts attribute that should be counted (default: transcript)
      --extraAttributes       FeatureCounts extra attriutes divided by comma (default: gene_name)
      --M                     FeatureCounts: Count multi-mapping reads
      --O                     FeatureCounts: Count reads overlapping features
      --fraction              FeatureCounts: Fractional counts for multi-mapping/overlapping features (must be used together with -M or -O)
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
params.noQualiRNA = false
params.noDepletionQC = false
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
genome_file = Channel.fromPath(params.reference, checkIfExists: true)
gff_file = Channel.fromPath(params.gff, checkIfExists: true)
rRNA_path = Channel.fromPath(params.rRNAgenes, checkIfExists:true)

isPaired = params.paired
noQualiControl = params.noQuali
noQualiControlRNA = params.noQualiRNA
noDepletionQC = params.noDepletionQC
noFC = params.noCounts
fcStrandness = params.featureCountsS
featureCounts_g = params.g
featureCounts_t = params.t
featureCounts_extra = params.extraAttributes
multiMapping = params.M
overlapping = params.O
fraction = params.fraction

// The basic input of the Pipeline. based on the Read name, a <base>_reference.info file is locatad, containing the information of the reference to be used for this sample
if(!isPaired){
    reads_ch = Channel.fromFilePairs(params.reads, size: 1, checkIfExists: true)
    reads_ch.ifEmpty{exit 1, 'Readfiles not specified'}
    reads_fastQC = reads_ch
    reads_trimgalore = reads_ch
} else{
    reads_ch = Channel.fromFilePairs(params.reads, size: 2, checkIfExists: true)
    reads_ch.ifEmpty{exit 1, "Readfiles not specified" }
    reads_fastQC = reads_ch
    reads_trimgalore = reads_ch
}

// Unzip genome file
process unZipGenome{
    conda 'environment.yml'
    input:
    file zipped_genome
    
    output:
    file "${zipped_genome.baseName}" 

    script:
    """
    gunzip -c $zipped_genome > ${zipped_genome.baseName}
    """
}
process unZipGFF{
    conda 'environment.yml'
    input:
    file zipped_gff
    
    output:
    file "${zipped_gff.baseName}" 

    script:
    """
    gunzip -c $zipped_gff > ${zipped_gff.baseName}
    """
}

// Convert GFF to GTF
process convertGFFtoGTF{
    conda 'environment.yml'
    tag "$gff"

    input:
    file gff

    output:
    file "${gff.baseName}.gtf"

    script:
    """
    gffread $gff -F -T -o ${gff.baseName}.gtf
    """
}

// FastQC for raw reads
process fastqc{
    conda 'environment.yml'
    tag "$id"
    label 'process_medium'

    publishDir "$pubDir/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    noQualiControl == false

    input: 
    tuple val(id) , file(reads)

    output:
    file '*_fastqc.{zip, html}'

    script:
    """
    fastqc --threads ${task.cpus} --quiet $reads
    """

}

// 1) Trimm the reads
process trimming{
    conda 'environment.yml'
    tag "$id"

    input:
    tuple val(id), file(read_file)

    publishDir "$pubDir/trimgalore/", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else null
        }
    
    output:
    tuple val(id), file("*.fq")
    file '*_fastqc.{zip,html}'
    file '*trimming_report.txt'

    script: 
    if(!isPaired){
        """
        trim_galore --dont_gzip --fastqc $read_file
        """
    }else{
        """
        trim_galore --dont_gzip --paired --fastqc ${read_file[0]} ${read_file[1]}
        """
    }
}

// 2) Create hisat_index
process hisat2_index{
    conda 'environment.yml'
    input:
    cpus 8
    file reference_genome

    output:
    file "${reference_genome.baseName}.*.ht2*"

    script:
    """
    hisat2-build -p ${task.cpus} $reference_genome ${reference_genome.baseName}.hisat2_index
    """
}

// 3) Create hisat_mapping
process hisat2_mapping{
    conda 'environment.yml'
    
    input:
    file indeces
    tuple val(id), file(reads)

    output:
    tuple val(id), file("*.bam")
    file("*bam")
    tuple val(id), file("*_summary.txt")

    script:
    index_name = indeces[0].toString() - ~/.\d.ht2l?/
    s = ""

    if(!isPaired){
        if(fcStrandness == "forward"){
          s = "F"
        } else if(fcStrandness == "reverse"){
          s = "R"
        }
      """
      hisat2 -x $index_name -U $reads --no-spliced-alignment --rna-strandness $s -p ${task.cpus} --met-stderr --new-summary --summary-file ${id}_summary.txt --rg-id ${id} --rg SM:${id} | samtools view -bS - > ${id}.bam
      """
    } else{
      """
      hisat2 -x $index_name --1 ${reads[0]} --2 ${reads[1]} -p ${task.cpus} --met-stderr --new-summary --summary-file ${id}_summary.txt --rg-id ${id} --rg SM:${id} | samtools view -bS - > ${id}.bam
      """
    }
}

// 4) counts features
process feature_counts{
    input:
    val g
    val t
    val extraAttributes
    file annotation
    file(bams)

    publishDir "$pubDir/FeatureCounts", mode: 'copy'

    when:
    noFC == false

    output:
    file("counts.txt")
    file("counts.txt.summary")

    script:
    M = ""
    O = ""
    frac = ""
    s = "0"
    p = ""
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
    if(isPaired){
      p ="-p"
    }
    """
    featureCounts -a $annotation -s $s -t $t -g $g --extraAttributes $extraAttributes -o counts.txt $M $O $frac $p $bams
    """
}

// samtools sort index
process samtools{
    tag "$id"
    cpus 4
    publishDir "$pubDir/Bams/", mode: 'copy'

    input:
    tuple val(id), file(bam) 

    output:
    tuple val(id), file("${id}.sorted.bam*")
    tuple val(id), file("${id}.sorted.bam*")
    

    script:
    """
    samtools sort -@ ${task.cpus} -o ${id}.sorted.bam $bam
    samtools index ${id}.sorted.bam
    """
}

//Qualimap on bam
process qualimap_bamqc{
    tag "$id"
    publishDir "$pubDir/QualiMapsBamQC/", mode: 'copy'

    when:
    noQualiControl == false && noFC == false && noQualiControlRNA == false

    input: 
    tuple val(id), file(bam)

    output:
    tuple val(id), file("${id}")

    script:
    sorted_bam = "${id}.sorted.bam"
    """
    qualimap bamqc -nt ${task.cpus} -bam $sorted_bam -outdir ${id}
    """
}

// Qualimap auf rnaseq
process qualimap_rnaseq{
    tag "$id"
    publishDir "$pubDir/QualiMapsRNAseq/", mode: 'copy'

    when:
    noQualiControl == false && noFC == false && noQualiControlRNA == false

    input:
    val annotation
    tuple val(id), file(bam)

    output:
    tuple val(id), file("${id}")

    script:
    sorted_bam = "${id}.sorted.bam"
    """
    qualimap rnaseq -bam $sorted_bam -gtf $annotation -outdir ${id}
    """
}

// MultiQC
process multiqc{
    conda 'environment.yml'
    publishDir "$pubDir/MultiQC", mode: 'copy'

    input:
    file(fastqc: 'fastqc/*')
    file('trimgalore/*')
    file('trimgaloreFastQC/*')
    file('alignment/*')
    file('qualimapBAMQC/*')
    file('qualimapRNAseq/*')
    file('featureCounts/*')

    output:
    file "*multiqc_report.html" 

    script:
    """
    multiqc .\\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m sortmerna -m fastqc -m qualimap -m salmon
    """
}

// Begin Depletion Analysis
// counts features for depletion analysis
process featureCountsAnalysis{
    publishDir "$pubDir/rRNADepletion", mode: 'copy'
    conda 'environment.yml'

    when:
    noDepletionQC == false

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
    script:
    """
    featureCounts -a $annotation -s $s -O -M --fraction -t $t -g $g --extraAttributes $extraAttributes -o counts.txt --minOverlap 10 --fracOverlap 0.9 $p $bams
    """
}

process depletionCalculationAnalysis{
    debug true
    conda 'environment.yml'
    publishDir "$pubDir/rRNADepletion", mode: 'copy'

    when:
    noDepletionQC == false

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

process genomecovAnalysis{
    conda 'environment.yml'
    publishDir "$pubDir/rRNADepletion", mode: 'copy'

    when:
    noDepletionQC == false

    input:
    each path(bams)

    output:
    path("*.bedgraph")

    script:
    """
    filename=\$(basename "$bams")
    filename_without_extension="\${filename%.*.*}"
    bedtools genomecov -d -ibam $bams > \$filename_without_extension.bedgraph
    """
 
}

process fasta_preparationAnalysis{
    conda 'environment.yml'
    publishDir "./bin", mode: 'copy'
    
    when:
    noDepletionQC == false

    input:
    path(fasta_file)

    output:
    path("fasta_data")

    script:
    """
    prepare_fasta_file.py -f $fasta_file
    """
}

process table_preparationAnalysis{
    conda 'environment.yml'
    publishDir "$pubDir/rRNADepletion", mode: 'copy'

    when:
    noDepletionQC == false

    input:
    path(genes_and_type)
    path(counts_file)
    path(gff_file)

    output:
    path('table_data.csv')

    script:
    """
    variance_calculation.py -fcf $counts_file -genesf $genes_and_type  -gff $gff_file
    """
}


workflow{
    genome_ch = unZipGenome(genome_file) //genome_fasta
    genome_analysis = genome_ch
    gff_ch = unZipGFF(gff_file)
    gff_depletion_analysis_1 = gff_ch
    gff_depletion_analysis_2 = gff_ch

    gtf_ch = convertGFFtoGTF(gff_ch)
    gtf_2 = gtf_ch
    gtf_depletion_analysis = gtf_ch
    fastqc_results_ch = fastqc(reads_fastQC) 
  
    trimming_ch = trimming(reads_trimgalore)
    trimmed_ch = trimming_ch[0]
    trimgalore_fastqc_results = trimming_ch[1]
    trimgalore_results = trimming_ch[2]

    hisat2_index_ch = hisat2_index(genome_ch)
    hisat2_mapping_ch = hisat2_mapping(hisat2_index_ch.collect(), trimmed_ch)
    alignment_files = hisat2_mapping_ch[0]
    alignment_files_fc = hisat2_mapping_ch[1]
    alignment_logs = hisat2_mapping_ch[2]

    featureCounts_ch = feature_counts(featureCounts_g, featureCounts_t, featureCounts_extra, gtf_ch, alignment_files_fc.collect())
    feture_counts = featureCounts_ch[0]
    featureCounts_logs = featureCounts_ch[1]

    samtools_ch = samtools(alignment_files)
    sorted_alignment_files_bamqc_ch = samtools_ch[0]
    alignment_files_analysis = samtools_ch[0]
    files_for_genomecov_analysis = samtools_ch[0]
    sorted_alignment_files_rnaseq_ch = samtools_ch[1]
    
    qualimap_bamqc_results_ch = qualimap_bamqc(sorted_alignment_files_bamqc_ch)
    qualimap_rnaseq_results_ch = qualimap_rnaseq(gtf_2.first(), sorted_alignment_files_rnaseq_ch)

    multiqc_report_ch = multiqc(fastqc_results_ch.collect().ifEmpty([]), 
    trimgalore_results.collect().ifEmpty([]),
    trimgalore_fastqc_results.collect().ifEmpty([]),
    alignment_logs.collect().ifEmpty([]),
    qualimap_bamqc_results_ch.collect().ifEmpty([]),
    qualimap_rnaseq_results_ch.collect().ifEmpty([]),
    featureCounts_logs.collect().ifEmpty([])
    )

    // Execute Analysis Pipeline
    
    rRNA_path_copy = rRNA_path
    bams = alignment_files_analysis.flatten().filter{file(it).isFile() && it.name.endsWith('.bam')}.collect()
    bams_genomecov = bams

    ch_feature_counts = featureCountsAnalysis(featureCounts_g, featureCounts_t, featureCounts_extra, gtf_depletion_analysis, bams)
    ch_feature_counts
        | flatten
        | branch {it ->
            txt: it.fileName.endsWith('counts.txt')
            summary:  it.fileName.endsWith('counts.txt.summary')
        }
        | set { result }

    counts_file = result.txt
    ch_depletion_calculation = depletionCalculationAnalysis(result.txt, featureCounts_g, rRNA_path)

    ch_genomecov = genomecovAnalysis(bams_genomecov)

    ch_fasta_file = fasta_preparationAnalysis(genome_analysis)

    ch_table = table_preparationAnalysis(rRNA_path_copy, counts_file, gff_depletion_analysis_2)
    

}