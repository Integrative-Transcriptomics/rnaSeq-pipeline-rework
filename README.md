# rnaSeq-Analysis

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --reads  "*.fastq.gz" --reference <.fna.gz> --gff <.gff.gz> -profile conda
    Mandatory arguments:
      --reads                 FastQ (.gz) input file containg the reads
      --reference             Input File in fasta format (.gz), nucleotide sequences.
      --gff                   Input file in gff (.gz) format, annotations
    Optional arguments:
      --paired                Paired end reads
      --noQuali               Disable quality control
      --noCounts              Disable feature counts
      --featureCountsS        FeatureCounts attribute strandedness: unstranded, reverse or forward strandness (default: reverse)
      --g                     FeatureCounts attribute to group features (default: locus_tag)
      --t                     FeatureCounts attribute that should be counted (default: transcript)
      --extraAttributes       FeatureCounts extra attriutes divided by comma (default: gene_name)
      --M                     FeatureCounts: Count multi-mapping reads
      --O                     FeatureCounts: Count reads overlapping features
      --fraction              FeatureCounts: Fractional counts for multi-mapping/overlapping features (must be used together with -M or -O)
      --pubDir                The directory where the results will be stored [def: Results]
