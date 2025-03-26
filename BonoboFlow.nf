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

def helpMessage() {
    log.info"""

To run BonoboFlow, use the following command:

 \
conda activate bonoboflow

nextflow run BonoboFlow.nf -resume  --in_fastq <directory to input files> \
                                    --outfile <directory to output files> \
                                    --ref_genome <directory to reference genome> \
                                    --sample_id <csv of sample IDs and barcode ID> \
                                    -w <directory to save the work files>

    

Mandatory arguments:
      --in_fastq                  Path to input fastq dirctory. Note: If you specify this you dont have to specify the  --raw_file
      --raw_file                  Path to raw POD5 or FAST5 files. Note: If you specify this, make sure you change the --basecalling flag to ON
      --outfile                   Path to output directory
      --ref_genome                reference sequence
      --sample_id                 a csv file containing barcode_Ids and sample Ids. An example csv file is provided in the BonoboFlow directory

Other arguments:
      --cpu                       cpus to used during the analyis. (default: 8)
      --phred                     minimum sequence quality score (default: 12)
      --lowerlength               set the lower length for input reads filter (default: 1000)
      --upperlength               set the upper length for input reads filter (default: 20000)
      --memory                    Memory allocated for each process. (default: 30 GB)
      --pipeline                  Specify whether you want to do genome assembly or haplotype reconstruction. (default: haplotype)
      --genomesize                Only required if you are running genome assembly (default: 5k)
      --basecalling               Please specify whether you would like to carry out basecalling (default: OFF). If "ON" ensure to provide raw files

Basecalling arguments:
      --basecallers               specify the basecalling tool (default: basecaller the alternative: duplex)
      --model                     Please specify the spped to run basecalling, the (default: sup the alternatives: fast, hac)

Barcoding arguments:
      --barcods                   barcods used during sequencing. (default: "EXP-NBD104 EXP-NBD114")
      --min_score_rear_barcode    minimum quality of of rear barcode (default: 75)
      --min_score_front_barcode   minimum quality of a rear barcode (default: 75)

mapping argements:
      --min_mq                    have mapping quality (default: 30)

Error_correction with rattle arguments:
      --repr_percentile           cluster representative percentile (default: 0.15)
      --score_threshold           minimum score for two reads to be in the same gene cluster (default: 0.2)
      --kmer_size                 k-mer size for isoform clustering (default: 11, maximum: 16)

Error_correction with vechat arguments:
      --split-size                split target sequences into chunks of desired size in lines (default: 5000)
      --cudapoa_batches           number of batches for CUDA accelerated polishing (default: 0)
      --cudaaligner-batches       number of batches for CUDA accelerated alignment (default: 0)
      

Haplotype arguments:
      --maxLD_floats              Maximum local divergence allowed for merging haplotypes. (default: 0.05)
      --maxGD_floats              Maximum global divergence allowed for merging haplotypes. (default: 0.05)
      --rmMisassembly_bool        Break contigs at potential misassembled positions (default: False)
      --correctErr_bool           Perform error correction for input reads (default: False)
      --minAbun_floats            Minimum abundance for filtering haplotypes (default: 0.02)
      --topks                     k seed reads size for haplotype construction (default: 100)
      --minovlplens               Minimum read overlap length. (default: 1000)
      --minseedlens               Minimum seed read length. (default: 2000)
      --maxohs                    Maximum overhang length allowed for overlaps. (default: 20)


""".stripIndent()
}

pipelineLogo()


// Show help emssage
params.help = false
if (params.help){
    helpMessage()
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
* Run dorado
*/

process runDorado {

    tag "raw_file"
    publishDir params.outfile, mode: "copy", overwrite: false
    
    input:
    path (raw_file)
    val (cpu)
    val (basecallers)
    val (model)
    
    output:
    path(basecalled)                
    path("basecalled"), emit: in_fastq
 
    script:
    """
    #!/usr/bin/env bash
    
    # Create the output directory
    mkdir -p basecalled
    
    # Run the dorado command
        dorado "$basecallers" --emit-fastq "$model" "$raw_file" > "basecalled/basecalled.fastq"
    done
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
    val (cpu)
    val (upperlength)
    val (phred)
    val (lowerlength)
    
    output:
    path choped_seq, emit: chopped
    
    script:
    """
    #!/usr/bin/env bash
    
    # Create the output directory
    mkdir -p choped_seq
    
    # Loop through input files and run Chopper
    for file in \$(ls $in_fastq); do
        sample_id=\$(basename "\$file")
        
        # Run the Chopper command
        chopper -q "$phred" -l "$lowerlength" --maxlength "$upperlength" --threads "$cpu" -i "$in_fastq/\$file" > "choped_seq/\$sample_id"
    done
    """
}

/*
* Run Barcoding
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
    path (demultiplexed)
    path (ref_genome)
    val (lowerlength)
    val (upperlength)
    val (cpu)

    output:
    path("mapped_reads"), emit: mapped

    script:
    """
    #!/usr/bin/env python

    import os
    import subprocess

    demultiplexed_dir = "${demultiplexed}"
    ref_gen = "${ref_genome}"  
    lowerlength = "${lowerlength}"
    upperlength = "${upperlength}"

    output_dir = os.path.join(os.path.dirname(demultiplexed_dir), "mapped_reads")
    os.makedirs(output_dir, exist_ok=True)

    # Copy the reference genome to the demultiplexed folder
    ref_genome_dest = os.path.join(demultiplexed_dir, os.path.basename(ref_gen))
    shutil.copyfile(ref_gen, ref_genome_dest)

    for dir in os.listdir(demultiplexed_dir):
        if not dir.startswith("bar"):
            continue

        barcode_id = dir
        dir_path = os.path.join(demultiplexed_dir, dir)

        if dir == 'barcoding_summary.txt':  # Skip the barcoding_summary.txt file
            continue

        output_subdir = os.path.join(output_dir, dir)
        os.makedirs(output_subdir, exist_ok=True)

        # Run the mapping 
        command = f"bwa index {ref_genome_dest} && cat {dir_path}/*.fastq | bwa mem -t ${cpu} {ref_genome_dest} /dev/stdin | samtools sort -@ ${cpu} | samtools view -F 4 -o {output_subdir}/{barcode_id}_align_sorted.bam"
        subprocess.run(command, shell=True, check=True)

        # Convert the BAM file to FASTQ
        bam_file = os.path.join(output_subdir, f"{barcode_id}_align_sorted.bam")
        fq_file = os.path.join(output_subdir, f"{barcode_id}_pre_mapped_reads.fastq")
        command = f"samtools bam2fq {bam_file} > {fq_file}"
        subprocess.run(command, shell=True, check=True)

        # Filter the reads based on length
        mapped_file = os.path.join(output_subdir, f"{barcode_id}_mapped_reads.fastq")
        command = f"filtlong --min_length {lowerlength} --max_length {upperlength} {fq_file} > {mapped_file}"
        subprocess.run(command, shell=True, check=True)
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

    mapped_dir = "${mapped}"
    output_dir = os.path.join(os.path.dirname(mapped_dir), "error_correction")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the mapped directory
    subdirs = [d for d in os.listdir(mapped_dir) if os.path.isdir(os.path.join(mapped_dir, d)) and d.startswith("bar")]

    for subdir in subdirs:
        subdir_path = os.path.join(mapped_dir, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Define the desired output file path
        corrected_output = os.path.join(output_subdir, f"{subdir}_corrected.fasta")

        # Run the error correction with direct output naming
        command1 = f"vechat {subdir_path}/{subdir}_mapped_reads.fastq --threads ${cpu} --platform ont --cudapoa-batches ${gpu} --cudaaligner-batches ${gpu} --split --split-size 5000 -o {corrected_output}"
        subprocess.run(command1, shell=True, check=True)
    """
}

/*
* Run genome_assembly
*/

process runAssembly {
    tag {"genome_assembly"}
    errorStrategy 'ignore'
    label 'bonobo_img'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path (corrected)
    val (genomesize)
    val (cpu)

    output:
    path("assembly"), emit: draftgenome

    script:
    """
    #!/usr/bin/env python

    import os
    import subprocess

    genomesize = "${genomesize}"
    corrected_dir = "${corrected}"
    output_dir = os.path.join(os.path.dirname(corrected_dir), "assembly")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the corrected directory
    subdirs = [d for d in os.listdir(corrected_dir) if os.path.isdir(os.path.join(corrected_dir, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(corrected_dir, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Identify the corrected read file
        corrected_read_file = None
        for filename in os.listdir(subdir_path):
            if filename.endswith("_corrected.fastq"):
                corrected_read_file = os.path.join(subdir_path, filename)
                converted_file = os.path.join(subdir_path, "reads.corrected.tmp2.fa")
                convert_cmd = f"seqtk seq -A {corrected_read_file} > {converted_file}"
                subprocess.run(convert_cmd, shell=True, check=True)
                corrected_read_file = converted_file
                break
            elif filename.endswith("_corrected.fasta"):
                corrected_read_file = os.path.join(subdir_path, filename)
                break

        if corrected_read_file is None:
            print(f"No corrected read file found in {subdir_path}. Skipping...")
            continue

        # Define output file paths
        assembled_file = os.path.join(output_subdir, f"{subdir}_contigs.fasta")
        final_assembly_file = os.path.join(output_subdir, "assembly.fasta")

        try:
            # Run the assembly
            command1 = f"flye --genome-size {genomesize} --threads {cpu} --nano-raw {corrected_read_file} --no-alt-contigs --scaffold --trestle --out-dir {output_subdir}"
            subprocess.run(command1, shell=True, check=True)

            # Check if assembly.fasta was created
            if os.path.exists(final_assembly_file):
                # If the file exists, move it
                command2 = f"mv {final_assembly_file} {assembled_file}"
                subprocess.run(command2, shell=True, check=True)
            else:
                # If the file does not exist, create a zero-byte file
                with open(assembled_file, 'w') as f:
                    pass

        except subprocess.CalledProcessError as e:
            print(f"Error occurred for {subdir}: {e}")
            # Create a zero-byte file in case of failure
            with open(assembled_file, 'w') as f:
                pass
    """
}


/*
* Run haplotype
*/
process runHaplotype {
    tag {"genome_haplotype"}
    errorStrategy 'ignore'
    label 'bonobo_img'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path (corrected)
    val (cpu)
    val (maxLD_floats)
    val (rmMisassembly_bool)
    val (correctErr_bool)
    val (minAbun_floats)
    val (maxGD_floats)
    val (topks)
    val (minovlplens)
    val (minseedlens)
    val (maxohs)

    output:
    path("haplotype"), emit: draftgenome

    script:
    """
    #!/usr/bin/env python

    import os
    import subprocess

    corrected_dir = "${corrected}"
    output_dir = os.path.join(os.path.dirname(corrected_dir), "haplotype")
    os.makedirs(output_dir, exist_ok=True)

    # Get the subdirectories in the corrected directory
    subdirs = [d for d in os.listdir(corrected_dir) if os.path.isdir(os.path.join(corrected_dir, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(corrected_dir, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Identify the corrected read file
        corrected_read_file = None
        for filename in os.listdir(subdir_path):
            if filename.endswith("_corrected.fastq"):
                corrected_read_file = os.path.join(subdir_path, filename)
                converted_file = os.path.join(subdir_path, "reads.corrected.tmp2.fa")
                convert_cmd = f"seqtk seq -A {corrected_read_file} > {converted_file}"
                subprocess.run(convert_cmd, shell=True, check=True)
                corrected_read_file = converted_file
                break
            elif filename.endswith("_corrected.fasta"):
                corrected_read_file = os.path.join(subdir_path, filename)
                break

        if corrected_read_file is None:
            print(f"No corrected read file found in {subdir_path}. Skipping...")
            continue

        # Define output file paths
        haplotype_file = os.path.join(output_subdir, f"{subdir}_contigs.fasta")
        final_haplotype_file = os.path.join(output_subdir, "haplotypes.final.fa")

        try:
            # Run the strainline command
            command1 = f"/app/Strainline/src/strainline.sh -i {corrected_read_file} -o {output_subdir} -p ont --maxLD ${maxLD_floats} --maxGD ${maxGD_floats} --rmMisassembly ${rmMisassembly_bool} --correctErr ${correctErr_bool} --minAbun ${minAbun_floats} -k ${topks} --minOvlpLen ${minovlplens} --minSeedLen ${minseedlens} --maxOH ${maxohs} --threads ${cpu}"
            subprocess.run(command1, shell=True, check=True)

            # Check if haplotypes.final.fa was created
            if os.path.exists(final_haplotype_file):
                # If the file exists, move it
                command2 = f"mv {final_haplotype_file} {haplotype_file}"
                subprocess.run(command2, shell=True, check=True)
            else:
                # If the file does not exist, create a zero-byte file
                with open(haplotype_file, 'w') as f:
                    pass

        except subprocess.CalledProcessError as e:
            print(f"Error occurred for {subdir}: {e}")
            # Create a zero-byte file in case of failure
            with open(haplotype_file, 'w') as f:
                pass

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
    path(pilon)
    
    output:
    path("proovframe"), emit: final_seq
   
    script: 
    """
    #!/usr/bin/env python

    import os
    import subprocess

    polished = "${pilon}"

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
                    {subdir_path}/{subdir}_polishing.fasta && \
                    proovframe fix -o {output_subdir}/{subdir}.fasta \
                    {subdir_path}/{subdir}_polishing.fasta \
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
    # Ensure the output directories exist
        mkdir -p mapped_reads/barcodex

    # Align reads with Minimap2, then sort with Samtools
        minimap2 -ax map-ont --eqx --MD --cs=long -t ${cpu} "${ref_genome}" ${chopped}/*.fastq | \
        samtools sort -@ ${cpu} -o mapped_reads/barcodex/barcodex_align_sorted.bam || { echo 'Minimap2 or Samtools sort failed'; exit 1; }

    # Extract IDs of mapped reads (without sequence modification)
        samtools view -F 4 mapped_reads/barcodex/barcodex_align_sorted.bam | cut -f1 | sort | uniq > mapped_reads/barcodex/mapped_read_ids.txt

    # Retrieve original reads using seqkit
        seqkit grep -f mapped_reads/barcodex/mapped_read_ids.txt ${chopped}/*.fastq > mapped_reads/barcodex/barcodex_pre_mapped_reads.fastq || { echo 'Seqkit filtering failed'; exit 1; }

    # Filter reads based on length using Filtlong
        filtlong --min_length ${lowerlength} --max_length ${upperlength} mapped_reads/barcodex/barcodex_pre_mapped_reads.fastq > mapped_reads/barcodex/barcodex_mapped_reads.fastq || { echo 'Filtlong filtering failed'; exit 1; }

    # Clean up intermediate files if needed
        rm mapped_reads/barcodex/barcodex_pre_mapped_reads.fastq
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
        if [[ -d "\$subdir" && "\$(basename \"\$subdir\")" == bar* ]]; then
            subdir_name="\$(basename \"\$subdir\")"
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
* Workflow
*/

workflow {
    if (params.basecalling == 'OFF') {
        runChopper(params.in_fastq, params.cpu, params.upperlength, params.phred, params.lowerlength)
    }

    else if (params.basecalling == 'ON') {
        runDorado(params.raw_file, params.cpu, params.basecallers, params.model)
        runChopper(runDorado.out.in_fastq, params.cpu, params.upperlength, params.phred, params.lowerlength)
        }

    if (params.demultiplexing == 'ON') {
        runBarcoding(runChopper.out.chopped, params.barcods, params.cpu, params.min_score_rear_barcode, params.min_score_front_barcode)
        runMapping(runBarcoding.out.demultiplexed, params.ref_genome, params.lowerlength, params.upperlength, params.cpu)
        
        // Run error correction based on selected tool
        if (params.error_correction_tool == 'vechat') {
            runErrcorrectVechat(runMapping.out.mapped, params.cpu, params.gpu)
        } else if (params.error_correction_tool == 'rattle') {
            runErrcorrectRattle(runMapping.out.mapped, params.cpu, params.repr_percentile, params.score_threshold, params.kmer_size)
        }

        // Get the appropriate corrected reads channel
        def corrected_reads = params.error_correction_tool == 'vechat' ? runErrcorrectVechat.out.corrected : runErrcorrectRattle.out.corrected
 
        if (params.pipeline == 'assembly') {
            runAssembly(corrected_reads, params.genomesize, params.cpu)
            runPolish_medaka(runAssembly.out.draftgenome, runMapping.out.mapped, params.cpu)
        }

         else if (params.pipeline == 'haplotype') {
            runHaplotype(corrected_reads, params.cpu, params.maxLD_floats, params.rmMisassembly_bool, params.correctErr_bool, params.minAbun_floats, params.maxGD_floats, params.topks, params.minovlplens, params.minseedlens, params.maxohs)
            runPolish_medaka(runHaplotype.out.draftgenome, runMapping.out.mapped, params.cpu)
        }
        runPolish_pilon(runPolish_medaka.out.medaka, runMapping.out.mapped, params.min_mq, params.cpu)
    }

    else if (params.demultiplexing == 'OFF') {
        runMapping_2(runChopper.out.chopped, params.ref_genome, params.lowerlength, params.upperlength, params.cpu)
        
        // Run error correction based on selected tool
        if (params.error_correction_tool == 'vechat') {
            runErrcorrectVechat(runMapping_2.out.mapped, params.cpu, params.gpu)
        } else if (params.error_correction_tool == 'rattle') {
            runErrcorrectRattle(runMapping_2.out.mapped, params.cpu, params.repr_percentile, params.score_threshold, params.kmer_size)
        }

        // Get the appropriate corrected reads channel
        def corrected_reads = params.error_correction_tool == 'vechat' ? runErrcorrectVechat.out.corrected : runErrcorrectRattle.out.corrected
    
        if (params.pipeline == 'assembly') {
            runAssembly(corrected_reads, params.genomesize, params.cpu)
            runPolish_medaka(runAssembly.out.draftgenome, runMapping_2.out.mapped, params.cpu)
        }

        else if (params.pipeline == 'haplotype') {
            runHaplotype(corrected_reads, params.cpu, params.maxLD_floats, params.rmMisassembly_bool, params.correctErr_bool, params.minAbun_floats, params.maxGD_floats, params.topks, params.minovlplens, params.minseedlens, params.maxohs)
            runPolish_medaka(runHaplotype.out.draftgenome, runMapping_2.out.mapped, params.cpu)
        }
        runPolish_pilon(runPolish_medaka.out.medaka, runMapping_2.out.mapped, params.min_mq, params.cpu)
    }

   runProovframe(runPolish_pilon.out.pilon)
   runSeqrenaming(runProovframe.out.final_seq, params.sample_id)
}