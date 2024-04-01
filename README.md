# INSTRUCTIONS FOR THE TOOL vis-rRNA-QC

Components:
Pipeline:       quality_control_and_depletion_analysis.nf
Visualization:  app.py



## Prerequisites
1)  Get the repository on your local machine
2)  Install the environment manager conda
    https://docs.conda.io/projects/conda/en/stable/
3)  Install Python 
4)  Install the folowing Python Packages with pip:
    - dash
    - dash-bootstrap-components
    - dash-core-components
    - dash-html-components
    - dash-table
    - pandas
    - plotly
    - numpy
5)  rRNA-genes file 'rRNA_genes_type.txt'
    Format:
        gene_id,type
        SCOr01,5S
        SCOr02,23S
        SCOr03,16S

## Pipeline Instructions
- Command for execution
  The command can only be executed if you are in the same folder as the pipeline quality_control_and_depletion_analysis.nf  
    **nextflow run quality_control_and_depletion_analysis.nf --reads "*fastq.gz" --reference <*.fna.gz> --gff <*.gff.gz> --rRNAgenes <rRNA_genes_type.txt> -profile conda**
- Additional parameters
    --paired                Paired end reads
    --noQuali               Disable quality control
    --noQualiRNA            Disable quality control RNASeq
    --noCounts              Disable feature counts
    --noDepletionQC         Disable depletion Analysis
    --featureCountsS        FeatureCounts attribute strandedness: unstranded, reverse or forward strandness (default: reverse)
    --g                     FeatureCounts attribute to group features (default: locus_tag)
    --t                     FeatureCounts attribute that should be counted (default: transcript)
    --extraAttributes       FeatureCounts extra attriutes divided by comma (default: gene_name)
    --M                     FeatureCounts: Count multi-mapping reads
    --O                     FeatureCounts: Count reads overlapping features
    --fraction              FeatureCounts: Fractional counts for multi-mapping/overlapping features (must be used
                            together with -M or -O)

## Visualization Instructions
- Command for execution:
  The command can only be executed if you are in the same folder as the app.py file, in this case you have
  to move to the 'bin' folder
    **python3 app.py**

  Output example:
  Dash is running on http://127.0.0.1:8050/

 * Serving Flask app 'app'
 * Debug mode: on

- After running th command above, the Dash app is running under the address displayed in the command line
  example:
    **http://127.0.0.1:8050/**
  Link can be opend with any browser
