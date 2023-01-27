#!/usr/bin/env nextflow
params.version = 1.0

def helpMessage() {
    log.info"""

Usage:
========================================================================================
                        B O N O B O F L O W     - P I P E L I N E
========================================================================================

  BonoboFlow PIPELINE for viral genome assembly pipeline from MinION sequenced reads

                ------------------------------------------------------

                      ------------------------------------------
                
                                --------------------
                
                                        -----
                            
                                          -

                      =============================================
                         BonoboFlow  ~  version ${params.version}
                      =============================================
    
                     To run BonoboFlow follow the command bellow:

nextflow run BonoboFlow.nf -resume --kit <sequencing kit> \
--flowcell <flow cell used during sequencing> \
--ref_genome <directory to reference genome> \
--input <directory to input files> \
--outfile <directory to output files> \
--mode <mode to run the pipeline> \
--processor <processor to be used during the analysis>
    

    Mandatory arguments:
      --kit                       Sequencit kit 
      --flowcell                  Flowcell used during sequencing
      --in_fast5                  Path to input fast5 dirctory 
      --outfile                   Path to output directory
      --ref_genome                reference sequence
      --mode                      mode for running the pipeline, defaults is FullPipeline, other mode include 
                                  OnlyGenomeAssembly (for doing only full/half/nearfull genome assemly)
                                  ORFprediction (for predicting coding region in the sequenced reads)
      --process                   process to be used.  the default is CPU

    Other arguments:
      --barcods                   barcods used during sequencing. The default barcoding kits are "EXP-NBD104 EXP-NBD114"
      --cpu                       cpus to used during the analyis. default is 8
      --memory                    memory allocated for each process. default in 30 GB
      --lowerlength               set the lower length for input reads filter (default: 150)
.
    """.stripIndent()
}


// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// validation of parameters 

ref = "${params.ref_genome}" 
File genome = new File(ref)

if (!genome.exists()){
    error "Error: Reference genome is not specified or you have provided incorrect path"
}

fast5 = "${params.in_fast5}" 
File input_seq = new File(fast5)

if (!input_seq.exists()){
    error "Error: Input fast5 path is not specified or you have provided incorrect path"
}


/*
* Run Basecalling
*/

process runBasecalling {
    tag {"fast5"}
    publishDir "$params.outfile", mode: "copy", overwrite: false
            
    input:
    path(in_fast5) 
    val(kit)
    val(flowcell)
    val(barcods)
    val(cpu)
    
    output:
    path(basecalled)
    path(demultiplexed), emit:  demultiplexed
    path("demultiplexed/*/*.fastq")
    path("basecalled/*.fastq")
    path("basecalled/*.log")
    path("basecalled/sequencing_summary.txt")


    script: 
    """
    mkdir basecalled demultiplexed
    guppy_basecaller -i ${in_fast5} \
        -s basecalled \
        --flowcell ${flowcell} \
        --kit ${kit} \
        --min_qscore 12 \
        --cpu_threads_per_caller "${params.cpu}"
    guppy_barcoder -i basecalled \
        -s demultiplexed \
        --trim_barcodes \
        --barcode_kits "${params.barcods}" \
        --require_barcodes_both_ends -t "${params.cpu}"
    """
}

/*
* Run Barcodng
*/
  
process runMapping {
    tag {"barcode_name"}
    label 'small_mem'
    publishDir "${params.outfile}", mode: "copy", overwrite: false
    
    input:
    path(demultiplexed) 
    val(ref_genome)
    val(lowerlength)
    
    output:
    path("demultiplexed/*/*mapped_reads.fastq"), emit:  mapped
    
    script: 
    """
    cd demultiplexed
    for dir in bar*/; do
		barcode_id=\${dir%*/}
		cat \${dir}/*.fastq | \
        minimap2 -a ${ref_genome} /dev/stdin | \
        samtools view -b | \
        samtools sort | \
        bamtools split -stub \${dir}/\${barcode_id}_align_sorted -mapped
        samtools bam2fq \${dir}/\${barcode_id}_align_sorted.MAPPED.bam > \
        \${dir}/\${barcode_id}_pre_mapped_reads.fastq
        filtlong --min_length ${lowerlength} \${dir}/\${barcode_id}_pre_mapped_reads.fastq > \
        \${dir}/\${barcode_id}_mapped_reads.fastq
    done
    """
}

/*
* Run error_correction
*/

process runErrcorrect {
    tag {"error_correction"}
    label 'small_mem'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(mapped)
    path(demultiplexed)
    
    output:
    path("demultiplexed/*/*_corrected.fastq"), emit:  correctedreads
   
    script: 
    """
    cd demultiplexed
    for dir in bar*/; do
		barcode_id=\${dir%*/}
        isONclust --consensus --t "${params.cpu}" --ont \
        --fastq \${dir}/\${barcode_id}_mapped_reads.fastq --outfolder \${dir}
        isONclust write_fastq --N 1 --clusters \${dir}/final_clusters.tsv \
        --fastq \${dir}/\${barcode_id}_mapped_reads.fastq --outfolder \${dir}/clusters
        run_isoncorrect --t 7  --fastq_folder \${dir}/clusters \
        --outfolder \${dir}/correction
        cat \${dir}/correction/*/*.fastq > \
        \${dir}/\${barcode_id}_corrected.fastq;
    done
    """
}

/*
* Run genome_assembly
*/

process runAssembly {
    tag {"genome_assembly"}
    publishDir "${params.outfile}", mode: "copy", overwrite: false


    input:
    path(demultiplexed)
    path(correctedreads)
    val(genomesize)
    
    output:
    path("demultiplexed/*/*unitigs.fasta"), emit:  draftgenome
   
    script: 
    """
    
    cd demultiplexed
    for dir in bar*/; do
		barcode_id=\${dir%*/}
        canu -d \${dir} \
             -p \${barcode_id} genomeSize=${genomesize} \
             -nanopore-corrected \${dir}/\${barcode_id}_corrected.fastq  \
             useGrid=false minReadLength=200  \
             minOverlapLength=100 \
             stopOnLowCoverage=0
    done
    """
}

/*
* Run polishing_medaka
*/

process runPolish_medaka {
    tag {"polish_medaka"}
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(demultiplexed)
    path(draftgenome)
    path(mapped)
    
    output:
    path("demultiplexed/*/*_consensus.fasta"), emit:  medaka_polishing
   
    script: 
    """
    cd demultiplexed
    for dir in bar*/; do
		barcode_id=\${dir%*/}
        medaka_consensus -i \${dir}/\${barcode_id}_mapped_reads.fastq \
                         -d \${dir}/\${barcode_id}.unitigs.fasta \
                         -o \${dir} \
                         -t 8
        mv \${dir}/consensus.fasta \
           \${dir}/\${barcode_id}_consensus.fasta
    done
    """
}

/*
* Run polishing_pilon
*/

process runPolish_pilon {
    tag {"polish_pilon"}
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(demultiplexed)
    path(medaka_polishing)
    path(mapped)
    
    output:
    path("demultiplexed/*/*_polishing.fasta"), emit:  pilon_polishing
    path("demultiplexed/*/*.vcf")
   
    script: 
    """
    cd demultiplexed
    for dir in bar*/; do
		barcode_id=\${dir%*/}  
        bwa index \${dir}/\${barcode_id}_consensus.fasta
        bwa mem -t7 \${dir}/\${barcode_id}_consensus.fasta \${dir}/\${barcode_id}_mapped_reads.fastq | \
        samtools sort  -@ 7 | \
        samtools markdup /dev/stdin \${dir}/\${barcode_id}_cons_ali_filter.bam
        samtools view -b -@ 7 -q 30 \${dir}/\${barcode_id}_cons_ali_filter.bam \
                      -o \${dir}/\${barcode_id}_cons_ali_dup_filter.bam
        samtools index -@ 7 \${dir}/\${barcode_id}_cons_ali_dup_filter.bam
        pilon -Xmx4096m -XX:-UseGCOverheadLimit \
              --genome \${dir}/\${barcode_id}_consensus.fasta \
              --unpaired \${dir}/\${barcode_id}_cons_ali_dup_filter.bam \
              --output \${dir}/\${barcode_id}_polishing \
              --vcf;
    done
    """
}

/*
* Run polishing_homopolish
*/

process runProovframe {
    label 'small_mem'
    tag {"proovframe"}
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(demultiplexed)
    path(pilon_polishing)
    
    output:
    path("demultiplexed/*/*_final_corrected.fasta"), emit: final_seq
    path(final_reports), emit: final_reports
   
   
    script: 
    """
    mkdir final_reports
    cd demultiplexed
    for dir in bar*/; do
		barcode_id=\${dir%*/}
        "${baseDir}"/packages/proovframe-main/bin/proovframe map \
        -a "${baseDir}"/db/viral_aa/viral_aa_ref_Seq.fasta \
        -o \${dir}/\${barcode_id}.tsv \
        \${dir}/\${barcode_id}_polishing.fasta
        "${baseDir}"/packages/proovframe-main/bin/proovframe fix \
        -o \${dir}/\${barcode_id}_final_corrected.fasta \
        \${dir}/\${barcode_id}_polishing.fasta \
        \${dir}/\${barcode_id}.tsv
    done
    """
}

/*
* Run renaming sequence
*/

process runSeqrenaming {
    label 'small_mem'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(demultiplexed)
    path(final_reports)
    path(final_seq)
    
    output:
    path("final_reports/*.fasta")   
   
    script: 
    """    
    cp demultiplexed/*/*_final_corrected.fasta final_reports
    """
}

workflow{
   runBasecalling (params.in_fast5, params.kit, params.flowcell, params.barcods, params.cpu)
   runMapping(runBasecalling.out.demultiplexed, params.ref_genome, params.lowerlength)   
   runErrcorrect(runBasecalling.out.demultiplexed, runMapping.out.mapped)
   runAssembly (runBasecalling.out.demultiplexed, runErrcorrect.out.correctedreads, params.genomesize)
   runPolish_medaka (runBasecalling.out.demultiplexed, runAssembly.out.draftgenome, runMapping.out.mapped)
   runPolish_pilon (runBasecalling.out.demultiplexed, runPolish_medaka.out.medaka_polishing, runMapping.out.mapped)
   runProovframe (runBasecalling.out.demultiplexed, runPolish_pilon.out.pilon_polishing)
   runSeqrenaming (runBasecalling.out.demultiplexed, runProovframe.out.final_seq, runProovframe.out.final_reports)
}