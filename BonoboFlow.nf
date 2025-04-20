#!/usr/bin/env nextflow

params.version = 1.0

def pipelineLogo() {
    log.info"""
========================================================================================
                        B O N O B O F L O W     - P I P E L I N E
========================================================================================

The BonoboFlow pipeline is a dedicated tool developed for the precise execution of viral 
       genome assembly and haplotypes construction from MinION sequencing reads

           ------------------------------------------------------------
 
                      ------------------------------------------

                               --------------------

                                      -----

                                        -

                   =============================================
                     BonoboFlow  ~  version ${params.version}
                   =============================================
""".stripIndent()
}

pipelineLogo()


// Show help emssage
params.help = false
if (params.help){
    log.info"""
BonoboFlow Pipeline v${params.version}
=====================================

The BonoboFlow pipeline is a dedicated tool for viral genome assembly and haplotype 
construction from MinION sequencing reads.

Usage:
    nextflow run BonoboFlow.nf [options]

Basic Command:
    nextflow run BonoboFlow.nf -resume \\
        --in_fastq <input_directory> \\
        --outfile <output_directory> \\
        --ref_genome <reference_genome> \\
        --sample_id <sample_csv_file> \\
        -w <work_directory>

Mandatory Arguments:
    --in_fastq                Path to input FASTQ directory (mutually exclusive with --raw_file)
    --raw_file               Path to raw POD5/FAST5 files (requires --basecalling ON)
    --outfile                Output directory path
    --ref_genome             Reference genome sequence
    --sample_id              CSV file with barcode and sample IDs

Optional Arguments:
    --cpu                    Number of CPUs to use (default: 8)
    --maxmem                 Memory allocation per process (default: 32 GB)
    --pipeline               Pipeline mode: "assembly" or "haplotype" (default: assembly)
    --phred                  Minimum sequence quality score (default: 12)
    --lowerlength            Minimum read length (default: 1000)
    --upperlength            Maximum read length (default: 20000)
    --genomesize             Expected genome size, assembly mode only (default: 5k)

Basecalling:
    --basecalling           Enable/disable basecalling (default: OFF)
    --basecallers           Tool choice: "basecaller" or "duplex" (default: basecaller)
    --model                 Model type: "sup", "fast", or "hac" (default: sup)

Barcoding:
    --barcods               Barcoding kits used (default: "EXP-NBD104 EXP-NBD114")
    --min_score_rear_barcode   Minimum rear barcode quality (default: 75)
    --min_score_front_barcode  Minimum front barcode quality (default: 75)

Error Correction:
    --error_correction_tool   Choose between "vechat" or "rattle" (default: vechat)
    --repr_percentile         Cluster representative percentile (default: 0.3)
    --score_threshold        Cluster similarity threshold (default: 0.2)
    --kmer_size             K-mer size for clustering (default: 12)

Haplotype Construction:
    --maxLD_floats          Maximum local divergence (default: 0.05)
    --maxGD_floats          Maximum global divergence (default: 0.05)
    --minAbun_floats        Minimum haplotype abundance (default: 0.2)
    --topks                 Seed reads size (default: 100)
    --minovlplens          Minimum overlap length (default: 1000)
    --minseedlens          Minimum seed length (default: 1000)
    --maxohs               Maximum overhang length (default: 20)

GPU Support:
    --gpu                   Enable GPU acceleration (default: 0)
                           Set to 1 for GPU support with compatible processes

For more information and examples, visit: https://github.com/nchis09/BonoboFlow
""".stripIndent()
    exit 0
}


// Parameter validation 

ref = "${params.ref_genome}" 
File genome = new File(ref)

if (!genome.exists()) {
    error "Error: Reference genome is not specified or you have provided an incorrect path"
}

sample_ids = "${params.sample_id}" 
File sample_id_csv = new File(sample_ids)

if (!sample_id_csv.exists()) {
    error "Error: Input sample ID path is not specified or you have provided an incorrect path"
}


def rawFile = "${params.raw_file}"
def inFastq = "${params.in_fastq}"

if ((rawFile && inFastq) || (!rawFile && !inFastq)) {
    error "Error: Only one of 'raw_file' or 'in_fastq' should be specified, not both"
}

if (rawFile) {
    def rawFileObj = new File(rawFile)
    if (!rawFileObj.exists()) {
        error "Error: Input raw file path is not specified or you have provided an incorrect path"
    }
} else {
    def inFastqObj = new File(inFastq)
    if (!inFastqObj.exists()) {
        error "Error: Input FASTQ file path is not specified or you have provided an incorrect path"
    }
}



/*
* Run dorado basecalling
*/

process runDoradobasecalling {
    tag "basecall"
    label 'dorado'
    publishDir params.outfile, mode: "copy", overwrite: false
    
    input:
    path (raw_file)
    val (cpu)
    val (basecallers)
    val (model)
    
    output:
    path "basecalled", emit: in_fastq
 
    script:
    """
    #!/usr/bin/env bash
    
    # Create the output directory
    mkdir -p basecalled
    
    # Run the dorado command
    dorado basecaller --device "cuda:all" --threads "${cpu}" "${model}" "${raw_file}" > "basecalled/basecalled.fastq"
    """
}

/*
* Run dorado demultiplexing
*/

process runDoradodemultiplexing {
    tag "demux"
    label 'dorado'
    publishDir params.outfile, mode: "copy", overwrite: false
    
    input:
    path (choped)
    val (cpu)
    val (barcods)
    
    output:
    path(demultiplexed_dir), emit: demultiplexed               

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Create the output directory
    mkdir -p demultiplexed_dir

    # Process both fastq and fq files
    echo "Processing input files..."
    
    # Try processing fastq files
    if ls "${choped}"/*.fastq 1> /dev/null 2>&1; then
        echo "Found .fastq files, processing..."
        dorado demux \\
            --threads "${cpu}" \\
            --emit-fastq \\
            --barcode-both-ends \\
            --kit-name "${barcods}" \\
            --output-dir "demultiplexed_dir" \\
            "${choped}"/*.fastq
    fi

    # Try processing fq files
    if ls "${choped}"/*.fq 1> /dev/null 2>&1; then
        echo "Found .fq files, processing..."
        dorado demux \\
            --threads "${cpu}" \\
            --emit-fastq \\
            --barcode-both-ends \\
            --kit-name "${barcods}" \\
            --output-dir "demultiplexed_dir" \\
            "${choped}"/*.fq
    fi

    # Organize fastq files by barcode, skip unclassified reads
    echo "Organizing demultiplexed files..."
    shopt -s nullglob  # Handle case when no files match pattern
    for fq in demultiplexed_dir/*.f*q; do
        [[ -f "\${fq}" ]] || continue  # Skip if not a file
        
        # Skip unclassified reads
        if [[ "\$(basename "\${fq}")" == *"unclassified"* ]]; then
            echo "Skipping unclassified reads: \${fq}"
            rm -f "\${fq}"
            continue
        fi

        # Extract barcode from filename
        barcode=\$(basename "\${fq}" | grep -o 'barcode[0-9]\\+' || echo "")
        if [[ -z "\${barcode}" ]]; then
            echo "Warning: Could not determine barcode for \${fq}"
            continue
        fi

        # Create barcode directory and move file
        mkdir -p "demultiplexed_dir/\${barcode}"
        mv "\${fq}" "demultiplexed_dir/\${barcode}/"
        echo "Moved \${fq} to demultiplexed_dir/\${barcode}/"
    done

    echo "Demultiplexing completed successfully"
    """
}

/*
* Run porewerchop
*/

process runPowerchop {

    tag "fastq"
    publishDir params.outfile, mode: "copy", overwrite: false
    label 'bonobo_img'
    
    input:
    path (in_fastq)
    val (cpu)
    
    output:
    path choped_seq, emit: choped
    
    script:
    """
    #!/usr/bin/env bash
    
    # Create the output directory
    mkdir -p choped_seq
    
    # Loop through input files and run porechop
    for file in \$(ls $in_fastq); do
        sample_id=\$(basename "\$file")
        
        # Run the porechop command
        porechop -i "$in_fastq/\$file" --threads "$cpu" -v 2 --extra_end_trim 0 --end_size 40 -o "choped_seq/\$sample_id"
    done
    """
}

/*
* Run Chopper
*/

process runChopper {
    tag "fastq"
    publishDir params.outfile, mode: "copy", overwrite: false
    label 'bonobo_img'
    
    input:
    path (in_fastq)
    val (phred)
    val (lowerlength)
    val (upperlength)
    val (cpu)

    output:
    path "choped_seq", emit: choped
    
    script:
    """
    #!/usr/bin/env bash
    
    # Create the output directory
    mkdir -p choped_seq
    
    # List all input files
    echo "Input directory contents:"
    ls -l "${in_fastq}"
    
    # Process all FASTQ files directly from the input directory
    for file in ${in_fastq}/*.fastq ${in_fastq}/*.fq; do
        # Check if file exists (to handle case when no .fq files exist)
        [ -e "\$file" ] || continue
        
        # Get just the filename without path
        filename=\$(basename "\$file")
        
        echo "Processing file: \$filename"
        
        # Run the Chopper command
        chopper -q "${phred}" -l "${lowerlength}" --maxlength "${upperlength}" --threads "${cpu}" -i "\$file" > "choped_seq/\$filename"
        
        # Check if chopper succeeded
        if [ \$? -ne 0 ]; then
            echo "Error: Chopper failed on file \$filename"
            exit 1
        fi
    done
    
    # Check if any files were processed
    shopt -s nullglob
    files=(choped_seq/*)
    if [ \${#files[@]} -eq 0 ]; then
        echo "No FASTQ files were processed in directory: ${in_fastq}"
        echo "Directory contents:"
        ls -la "${in_fastq}/"
        exit 1
    fi
    
    echo "Successfully processed files:"
    ls -l choped_seq/
    """
}

/*
* Run Barcoding with guppy
*/
process runBarcoding {
    tag "barcodes"
    publishDir params.outfile, mode: "copy", overwrite: false
            
    input:
    path (choped) 
    val (barcods)
    val (cpu)
    val (min_score_rear_barcode)
    val (min_score_front_barcode)
    
    output:
    path(demultiplexed_dir), emit: demultiplexed
    path("demultiplexed_dir/*/*.fastq")

    script: 
    """
    mkdir demultiplexed_dir
    guppy_barcoder -i "${choped}" -s demultiplexed_dir \
    --barcode_kits "${barcods}" -t "${cpu}" --fastq_out --min_score_rear_override "${min_score_rear_barcode}" --min_score "${min_score_front_barcode}" --trim_barcodes
    find demultiplexed_dir -maxdepth 1 -type d -exec du -s {} \\; | awk '\$1 < 1000 {print \$2}' | \
    xargs -I {} sh -c 'rm -rf {} && echo {} deleted' > demultiplexed_dir/deleted_directories.txt
    """
}

/*
* Run Mapping
*/

process runMapping {
    tag { "barcode_name" }
    label 'bonobo_img'
    publishDir path: "${params.outfile}", mode: "copy", overwrite: false

    input:
    path demultiplexed
    path ref_genome
    val lowerlength
    val upperlength
    val cpu

    output:
    path "mapped_reads", emit: mapped

    script:
    """
    #!/usr/bin/env bash
    set -e

    # Create output directory
    mkdir -p mapped_reads

    # Copy reference genome to working directory
    cp "${ref_genome}" ./reference.fasta

    for dir in ${demultiplexed}/bar*/; do
        if [ ! -d "\$dir" ]; then
            echo "No barcode directories found"
            exit 0
        fi
        
        barcode_id=\$(basename \${dir%/})
        echo "Processing barcode: \$barcode_id"
        
        # Check if directory contains fastq files
        if ! ls \${dir}/*.fastq >/dev/null 2>&1; then
            echo "No FASTQ files found in \$dir, skipping..."
            continue
        fi
        
        mkdir -p "mapped_reads/\${barcode_id}"
        
        # Concatenate all fastq files, map with minimap2, convert to bam, sort, and split
        echo "Mapping reads for \$barcode_id..."
        cat \${dir}/*.fastq | \
        minimap2 -a -t ${cpu} ./reference.fasta /dev/stdin | \
        samtools view -b | \
        samtools sort -@ ${cpu} | \
        bamtools split -stub "mapped_reads/\${barcode_id}/\${barcode_id}_align_sorted" -mapped

        # Check if the mapped BAM file exists and is not empty
        mapped_bam="mapped_reads/\${barcode_id}/\${barcode_id}_align_sorted.MAPPED.bam"
        if [ ! -s "\$mapped_bam" ]; then
            echo "No mapped reads found for \$barcode_id, removing directory..."
            rm -rf "mapped_reads/\${barcode_id}"
            continue
        fi
        
        # Convert mapped BAM to FASTQ
        echo "Converting BAM to FASTQ for \$barcode_id..."
        samtools bam2fq "\$mapped_bam" > \
        "mapped_reads/\${barcode_id}/\${barcode_id}_pre_mapped_reads.fastq"

        # Check if the converted FASTQ is empty
        if [ ! -s "mapped_reads/\${barcode_id}/\${barcode_id}_pre_mapped_reads.fastq" ]; then
            echo "No reads in FASTQ for \$barcode_id, removing directory..."
            rm -rf "mapped_reads/\${barcode_id}"
            continue
        fi
        
        # Filter reads by length
        echo "Filtering reads by length for \$barcode_id..."
        filtlong --min_length ${lowerlength} --max_length ${upperlength} "mapped_reads/\${barcode_id}/\${barcode_id}_pre_mapped_reads.fastq" > \
        "mapped_reads/\${barcode_id}/\${barcode_id}_mapped_reads.fastq"

        # Check if the final filtered FASTQ is empty
        if [ ! -s "mapped_reads/\${barcode_id}/\${barcode_id}_mapped_reads.fastq" ]; then
            echo "No reads passed length filtering for \$barcode_id, removing directory..."
            rm -rf "mapped_reads/\${barcode_id}"
            continue
        fi
        
        # Cleanup intermediate files
        rm -f "mapped_reads/\${barcode_id}/\${barcode_id}_pre_mapped_reads.fastq"
        
        echo "Completed processing for \$barcode_id"
    done
    
    # Clean up reference copy
    rm -f reference.fasta
    
    echo "All barcodes processed successfully"
    """
}

/*
* Run error_correction_vechat
*/

process runErrcorrectVechat {
    tag "error_correction"
    errorStrategy 'ignore'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path (mapped)
    val (cpu)
    val (gpu)

    output:
    path("error_correction"), emit: corrected

    script:
    """
    #!/usr/bin/env python

    import os
    import subprocess
    import glob
    import shutil

    mapped_dir = "${mapped}"
    output_dir = os.path.join(os.path.dirname(mapped_dir), "error_correction")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the mapped directory
    subdirs = [d for d in os.listdir(mapped_dir) if os.path.isdir(os.path.join(mapped_dir, d)) and d.startswith("bar")]

    for subdir in subdirs:
        subdir_path = os.path.join(mapped_dir, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Find the *_mapped_reads.fastq file
        mapped_files = glob.glob(os.path.join(subdir_path, "*_mapped_reads.fastq"))
        if not mapped_files:
            print(f"No *_mapped_reads.fastq file found in {subdir_path}. Skipping...")
            continue
        
        mapped_read_file = mapped_files[0]  # Use the first matching file

        # Run the error correction with direct output naming
        command1 = f"vechat {mapped_read_file} --threads ${cpu} --platform ont --cudapoa-batches ${gpu} --cudaaligner-batches ${gpu} --split --split-size 1000 -o {output_subdir}"
        subprocess.run(command1, shell=True, check=True)

        # Move and rename the corrected output file
        tmp_file = os.path.join(output_subdir, "reads.corrected.tmp2.fa")
        if os.path.exists(tmp_file):
            final_file = os.path.join(output_subdir, f"{subdir}_corrected.fastq")
            shutil.move(tmp_file, final_file)
            print(f"Moved corrected reads to {final_file}")
        else:
            print(f"Warning: Expected output file {tmp_file} not found")
    """
}

/*
* Run polishing_medaka
*/

process runPolish_medaka {
    tag {"polish_medaka"}
    errorStrategy 'ignore'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path (draftgenome)
    path (mapped)
    val (cpu)
    
    output:
    path("medaka"), emit: medaka
   
    script: 
    """
    #!/usr/bin/env python

    import os
    import subprocess

    mapped_dir = "${mapped}"
    assembled_dir = "${draftgenome}"

    output_dir = os.path.join(os.path.dirname(assembled_dir), "medaka")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the assembled directory
    subdirs = [d for d in os.listdir(assembled_dir) if os.path.isdir(os.path.join(assembled_dir, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(assembled_dir, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Get the corresponding mapped file
        mapped_file = os.path.join(mapped_dir, subdir, f"{subdir}_mapped_reads.fastq")

        # Define output file paths
        polish_medaka_file = os.path.join(output_subdir, f"{subdir}_polish_medaka.fasta")
        consensus_file = os.path.join(output_subdir, "consensus.fasta")

        try:
            # Run the medaka consensus command
            command1 = f"medaka_consensus -i {mapped_file} \
                             -d {subdir_path}/{subdir}_contigs.fasta \
                             -o {output_subdir} \
                             -m r10_min_high_g340 \
                             -t ${cpu}"
            subprocess.run(command1, shell=True, check=True)

            # Check if the consensus.fasta file was created
            if os.path.exists(consensus_file):
                # If the file exists, move it to the expected output file
                command2 = f"mv {consensus_file} {polish_medaka_file}"
                subprocess.run(command2, shell=True, check=True)
            else:
                # If the file does not exist, create a zero-byte placeholder
                with open(polish_medaka_file, 'w') as f:
                    pass

        except subprocess.CalledProcessError as e:
            print(f"Error occurred for {subdir}: {e}")
            # Create a zero-byte file in case of failure
            with open(polish_medaka_file, 'w') as f:
                pass

    """
}


/*
* Run polishing_pilon
*/

process runPolish_pilon {
    tag {"polish_pilon"}
    errorStrategy 'ignore'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path (medaka)
    path (mapped)
    val (cpu)
    val (min_mq)
    
    output:
    path("pilon"), emit: pilon
   
    script: 
    """
    #!/usr/bin/env python

    import os
    import subprocess

    mapped_dir = "${mapped}"
    polished = "${medaka}"

    output_dir = os.path.join(os.path.dirname(polished), "pilon")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the polished directory
    subdirs = [d for d in os.listdir(polished) if os.path.isdir(os.path.join(polished, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(polished, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Get the corresponding mapped file
        mapped_file = os.path.join(mapped_dir, subdir, f"{subdir}_mapped_reads.fastq")

        # Define the output pilon files
        pilon_output_prefix = os.path.join(output_subdir, f"{subdir}_polishing")
        pilon_vcf_file = f"{pilon_output_prefix}.vcf"

        try:
            # Run the pilon polishing pipeline
            command1 = f"bwa index {subdir_path}/{subdir}_polish_medaka.fasta && \
                         bwa mem -t ${cpu} {subdir_path}/{subdir}_polish_medaka.fasta {mapped_file} | \
                         samtools sort -@ ${cpu} | \
                         samtools markdup /dev/stdin {output_subdir}/{subdir}_cons_ali_filter.bam && \
                         samtools view -b -@ ${cpu} -q ${min_mq} {output_subdir}/{subdir}_cons_ali_filter.bam \
                                              -o {output_subdir}/{subdir}_cons_ali_dup_filter.bam && \
                         samtools index -@ ${cpu} {output_subdir}/{subdir}_cons_ali_dup_filter.bam && \
                         pilon -Xmx4096m -XX:-UseGCOverheadLimit \
                               --genome {subdir_path}/{subdir}_polish_medaka.fasta \
                               --unpaired {output_subdir}/{subdir}_cons_ali_dup_filter.bam \
                               --output {pilon_output_prefix} \
                               --vcf"
            subprocess.run(command1, shell=True, check=True)

            # Check if the pilon VCF file was created
            if os.path.exists(pilon_vcf_file):
                print(f"Pilon polishing completed successfully for {subdir}.")
            else:
                # If the VCF file does not exist, create a zero-byte placeholder
                with open(pilon_vcf_file, 'w') as f:
                    pass

        except subprocess.CalledProcessError as e:
            print(f"Error occurred for {subdir}: {e}")
            # Create a zero-byte placeholder file in case of failure
            with open(pilon_vcf_file, 'w') as f:
                pass
    """
}


/*
* Run polishing_prooveframe
*/


process runProovframe {
    label 'bonobo_img'
    tag {"proovframe"}
    errorStrategy 'ignore'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(medaka)
    
    output:
    path("proovframe"), emit: final_seq
   
    script: 
    """
    #!/usr/bin/env python

    import os
    import subprocess

    polished = "${medaka}"

    output_dir = os.path.join(os.path.dirname(polished), "proovframe")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the polished directory
    subdirs = [d for d in os.listdir(polished) if os.path.isdir(os.path.join(polished, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(polished, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Run the proovframe 
        command1 = f"proovframe map -a /app/db/viral_aa/viral_aa_ref_Seq.fasta \
                    -o {output_subdir}/{subdir}.tsv \
                    {subdir_path}/{subdir}_polish_medaka.fasta && \
                    proovframe fix -o {output_subdir}/{subdir}.fasta \
                    {subdir_path}/{subdir}_polish_medaka.fasta \
                    {output_subdir}/{subdir}.tsv"
        subprocess.run(command1, shell=True, check=True)
    """
}


/*
* Run renaming sequence
*/

process runSeqrenaming {
    label 'small_mem'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(proovframe)
    path(sample_id)
    
    output:
    path("final_reports/*.fasta")
   
    script:
    """
    #!/usr/bin/env python

    import os
    import csv
    import shutil

    # Create the final_reports directory
    os.makedirs("final_reports", exist_ok=True)

    # Get the subdirectories within the proovframe directory
    subdirs = [d for d in os.listdir("${proovframe}") if os.path.isdir(os.path.join("${proovframe}", d))]

    # Copy the fasta files to the final_reports directory
    for subdir in subdirs:
        subdir_path = os.path.join("${proovframe}", subdir)
        fasta_files = [f for f in os.listdir(subdir_path) if os.path.isfile(os.path.join(subdir_path, f)) and f.endswith('.fasta')]

        for fasta_file in fasta_files:
            fasta_file_path = os.path.join(subdir_path, fasta_file)
            new_file_path = os.path.join(os.getcwd(), "final_reports", fasta_file)
            shutil.copy(fasta_file_path, new_file_path)

    # Read the CSV file
    mapping = {}
    with open("${sample_id}", "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            barcodeID = row[0].strip()
            newSampleID = row[1].strip()
            mapping[barcodeID] = newSampleID

    # Rename the copied FASTA files
    final_reports_dir = os.path.join(os.getcwd(), "final_reports")
    fasta_files = [f for f in os.listdir(final_reports_dir) if os.path.isfile(os.path.join(final_reports_dir, f)) and f.endswith('.fasta')]

    for fasta_file in fasta_files:
        barcodeID = fasta_file.split(".")[0]
        newSampleID = mapping.get(barcodeID)

        if newSampleID is not None:
            new_file_name = f"{newSampleID}.fasta"
            new_file_path = os.path.join(final_reports_dir, new_file_name)

            # Rename the file
            os.rename(os.path.join(final_reports_dir, fasta_file), new_file_path)
    """
}

/*
* Run Mapping withought demultiplexing
*/

process runMapping_2 {
    tag 'mapping'
    label 'bonobo_img'
    publishDir path: "${params.outfile}", mode: "copy", overwrite: false

    input:
    path chopped
    path ref_genome
    val lowerlength
    val upperlength
    val cpu

    output:
    path "mapped_reads", emit: mapped

    script:
    """
    #!/usr/bin/env bash
    set -e

    # Create output directory
    mkdir -p mapped_reads/barcodex

    # Copy reference genome to working directory
    cp "${ref_genome}" ./reference.fasta

    # Concatenate all fastq files, map with minimap2, convert to bam, sort, and split
    cat ${chopped}/*.fastq | \
    minimap2 -a -t ${cpu} ./reference.fasta /dev/stdin | \
    samtools view -b | \
    samtools sort -@ ${cpu} | \
    bamtools split -stub mapped_reads/barcodex/barcodex_align_sorted -mapped
    
    # Convert mapped BAM to FASTQ
    samtools bam2fq mapped_reads/barcodex/barcodex_align_sorted.MAPPED.bam > \
    mapped_reads/barcodex/barcodex_pre_mapped_reads.fastq
    
    # Filter reads by length
    filtlong --min_length ${lowerlength} --max_length ${upperlength} mapped_reads/barcodex/barcodex_pre_mapped_reads.fastq > \
    mapped_reads/barcodex/barcodex_mapped_reads.fastq
    
    # Check if the final filtered FASTQ is empty
    if [ ! -s "mapped_reads/barcodex/barcodex_mapped_reads.fastq" ]; then
        echo "No reads passed length filtering, removing directory..."
        rm -rf "mapped_reads/barcodex"
    fi
    
    # Clean up reference copy
    rm -f reference.fasta
    """
}


/*
* Run error_correction_rattle
*/

process runErrcorrectRattle {
    tag "error_correction"
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path mapped
    val cpu
    val repr_percentile
    val score_threshold
    val kmer_size

    output:
    path "error_correction", emit: corrected

    script:
    """
    mkdir -p error_correction

    for subdir in $mapped/*; do
        if [[ -d "\$subdir" && "\$(basename "\$subdir")" == bar* ]]; then
            subdir_name="\$(basename "\$subdir")"
            output_subdir="error_correction/\${subdir_name}"
            mkdir -p "\${output_subdir}"

            rattle cluster -i "\${subdir}/\${subdir_name}_mapped_reads.fastq" \\
                -p $repr_percentile \\
                -s $score_threshold \\
                -t $cpu \\
                --iso-kmer-size $kmer_size \\
                -o "\${output_subdir}"

            rattle correct -i "\${subdir}/\${subdir_name}_mapped_reads.fastq" \\
                -c "\${output_subdir}/clusters.out" \\
                -t $cpu \\
                -o "\${output_subdir}"

            mv "\${output_subdir}/corrected.fq" "\${output_subdir}/\${subdir_name}_corrected.fastq"
        fi
    done
    """
}


/*
* Run genome_assembly
*/

process runAssembly {

    tag { "genome_assembly" }
    errorStrategy 'ignore'
    label 'bonobo_img'
    publishDir "${params.outfile ?: 'results'}", mode: "copy", overwrite: false

    input:
    path corrected
    val cpu 
    val genomesize 

    output:
    path("${params.output_dir ?: 'assembly'}"), emit: draftgenome

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    output_dir="${params.output_dir ?: 'assembly'}"
    mkdir -p "\${output_dir}"

    echo "Processing directory: ${corrected}"
    echo "Output directory: \${output_dir}"

    run_command() {
        local cmd="\$1"
        local desc="\$2"
        echo "Running: \${desc:-\$cmd}"
        if eval "\${cmd}"; then
            echo "Command succeeded."
        else
            echo "Error running: \${desc:-\$cmd}" >&2
            return 1
        fi
    }

    for subdir in "${corrected}"/bar*/; do
        [[ -d "\$subdir" ]] || continue

        barcode_name=\$(basename "\${subdir}")
        echo -e "\\nProcessing barcode: \${barcode_name}"

        corrected_read_file=\$(find "\${subdir}" -maxdepth 1 -name "*_corrected.fastq" | head -n 1)

        if [[ -z "\${corrected_read_file}" ]]; then
            echo "No corrected reads found in \${barcode_name}, skipping..."
            continue
        fi

        echo "Found corrected reads: \${corrected_read_file}"

        output_subdir="\${output_dir}/\${barcode_name}"
        mkdir -p "\${output_subdir}"

        fasta_file="\${output_subdir}/\${barcode_name}_corrected.fasta"
        if ! run_command "seqtk seq -A \${corrected_read_file} > \${fasta_file}" "Converting FASTQ to FASTA"; then
            echo "FASTQ to FASTA conversion failed for \${barcode_name}, skipping..."
            continue
        fi

        flye_cmd="flye --genome-size ${genomesize} \\
            --threads ${cpu} \\
            --nano-corr \${corrected_read_file} \\
            --no-alt-contigs \\
            --scaffold \\
            --trestle \\
            --out-dir \${output_subdir}"

        if run_command "\${flye_cmd}" "Running Flye"; then
            final_file="\${output_subdir}/assembly.fasta"
            output_file="\${output_subdir}/\${barcode_name}_contigs.fasta"
            if [[ -f "\${final_file}" ]]; then
                echo "Moving output to: \${output_file}"
                mv "\${final_file}" "\${output_file}"
            else
                echo "Warning: Flye output not found at \${final_file}"
                touch "\${output_file}"
            fi
        else
            echo "Flye failed for \${barcode_name}, creating empty output"
            touch "\${output_subdir}/\${barcode_name}_contigs.fasta"
        fi
    done

    echo -e "\\nGenome assembly process completed"

    """
}


/*
* Run haplotype
*/
process runHaplotype {
    tag { "genome_haplotype" }
    label 'bonobo_img'
    errorStrategy 'ignore'
    time '24h'       // Allow up to 24 hours runtime
    publishDir "${params.outfile}", mode: 'copy', overwrite: false

    input:
    path mapped
    val maxLD_floats
    val maxGD_floats
    val rmMisassembly_bool
    val minAbun_floats
    val topks
    val minovlplens
    val minseedlens
    val maxohs
    val cpu
    val maxmem

    output:
    path "haplotype", emit: draftgenome

    script:
    """
    mkdir -p haplotype

    for subdir in "${mapped}"/*; do
        if [[ ! -d "\$subdir" || \$(basename "\$subdir") != bar* ]]; then
            continue
        fi

        subname=\$(basename "\$subdir")
        haplo_out="haplotype/\$subname"
        mkdir -p "\$haplo_out"

        mapped_read=""
        for file in "\$subdir"/*; do
            if [[ "\$file" == *_mapped_reads.fastq ]]; then
                fasta_out="\$subdir/\${subname}_mapped_reads.fasta"
                if command -v seqtk >/dev/null 2>&1; then
                    seqtk seq -A "\$file" > "\$fasta_out"
                    mapped_read="\$fasta_out"
                else
                    echo "seqtk not found. Cannot convert fastq to fasta. Skipping \$subdir..." >&2
                    continue 2
                fi
                break
            elif [[ "\$file" == *_corrected.fasta ]]; then
                mapped_read="\$file"
                break
            fi
        done

        if [[ -z "\$mapped_read" ]]; then
            echo "No mapped read file found in \$subdir. Skipping..." >&2
            continue
        fi

        if ! timeout 1h /app/Strainline/src/strainline.sh \\
            -i "\$mapped_read" \\
            -p ont \\
            --maxLD ${maxLD_floats} \\
            --maxGD ${maxGD_floats} \\
            --rmMisassembly ${rmMisassembly_bool} \\
            --correctErr True \\
            --minAbun ${minAbun_floats} \\
            -k ${topks} \\
            --minOvlpLen ${minovlplens} \\
            --minSeedLen ${minseedlens} \\
            --maxOH ${maxohs} \\
            --threads ${cpu} \\
            --maxMem ${maxmem} \\
            -o "\$haplo_out" 2>&1; then
            status=\$?
            if [ \$status -eq 124 ]; then
                echo "Strainline timed out after 1 hour for \$subname" >&2
            else
                echo "Strainline failed with status \$status for \$subname" >&2
            fi
            continue
        fi

        if [[ -f "\$haplo_out/haplotypes.final.fa" ]]; then
            mv "\$haplo_out/haplotypes.final.fa" "\$haplo_out/\${subname}_contigs.fasta"
        else
            echo "No haplotype output for \$subname. Creating empty contigs file." >&2
            touch "\$haplo_out/\${subname}_contigs.fasta"
        fi
    done
    """
}

/*
* Workflow
*/

workflow {
    if (params.basecalling == 'OFF') {
        runChopper(params.in_fastq, params.phred, params.lowerlength, params.upperlength, params.cpu)
    }

    else if (params.basecalling == 'ON') {
        runDoradobasecalling(params.raw_file, params.cpu, params.basecallers, params.model)
        runChopper(runDoradobasecalling.out.in_fastq, params.phred, params.lowerlength, params.upperlength, params.cpu)
        }

    if (params.demultiplexing == 'ON') {
        runDoradodemultiplexing(runChopper.out.choped, params.cpu, params.barcods)
        runMapping(runDoradodemultiplexing.out.demultiplexed, params.ref_genome, params.lowerlength, params.upperlength, params.cpu)
        
        if (params.pipeline == 'assembly') {
            // Run error correction based on selected tool
            if (params.error_correction_tool == 'vechat') {
                runErrcorrectVechat(runMapping.out.mapped, params.cpu, params.gpu)
            } else if (params.error_correction_tool == 'rattle') {
                runErrcorrectRattle(runMapping.out.mapped, params.cpu, params.repr_percentile, params.score_threshold, params.kmer_size)
            }

            // Get the appropriate corrected reads channel
            def corrected_reads = params.error_correction_tool == 'vechat' ? runErrcorrectVechat.out.corrected : runErrcorrectRattle.out.corrected
            runAssembly(corrected_reads, params.cpu, params.genomesize)
            runPolish_medaka(runAssembly.out.draftgenome, runMapping.out.mapped, params.cpu)
            runProovframe(runPolish_medaka.out.medaka)
        }
        else if (params.pipeline == 'haplotype') {
            // Skip error correction for haplotype pipeline
            runHaplotype(runMapping.out.mapped, params.maxLD_floats, params.maxGD_floats, params.rmMisassembly_bool, params.minAbun_floats, params.topks, params.minovlplens, params.minseedlens, params.maxohs, params.cpu, params.maxmem)
            runPolish_medaka(runHaplotype.out.draftgenome, runMapping.out.mapped, params.cpu)
            runProovframe(runPolish_medaka.out.medaka)
        }
    }

    else if (params.demultiplexing == 'OFF') {
        runMapping_2(runChopper.out.choped, params.ref_genome, params.lowerlength, params.upperlength, params.cpu)
        
        if (params.pipeline == 'assembly') {
            // Run error correction based on selected tool
            if (params.error_correction_tool == 'vechat') {
                runErrcorrectVechat(runMapping_2.out.mapped, params.cpu, params.gpu)
            } else if (params.error_correction_tool == 'rattle') {
                runErrcorrectRattle(runMapping_2.out.mapped, params.cpu, params.repr_percentile, params.score_threshold, params.kmer_size)
            }

            // Get the appropriate corrected reads channel
            def corrected_reads = params.error_correction_tool == 'vechat' ? runErrcorrectVechat.out.corrected : runErrcorrectRattle.out.corrected
            runAssembly(corrected_reads, params.cpu, params.genomesize)
            runPolish_medaka(runAssembly.out.draftgenome, runMapping_2.out.mapped, params.cpu)
            runProovframe(runPolish_medaka.out.medaka)
        }
        else if (params.pipeline == 'haplotype') {
            // Skip error correction for haplotype pipeline
            runHaplotype(runMapping_2.out.mapped, params.maxLD_floats, params.maxGD_floats, params.rmMisassembly_bool, params.minAbun_floats, params.topks, params.minovlplens, params.minseedlens, params.maxohs, params.cpu, params.maxmem)
            runPolish_medaka(runHaplotype.out.draftgenome, runMapping_2.out.mapped, params.cpu)
            runProovframe(runPolish_medaka.out.medaka)
        }
    }
   runSeqrenaming(runProovframe.out.final_seq, params.sample_id)
}

workflow.onComplete {
    println "="*80
    println "Final Consensus Sequences Are Saved In:"
    println "${params.outfile}/final_reports"
    println "="*80
}