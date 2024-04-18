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
conda activate nextflow

nextflow run BonoboFlow.nf -resume \
--ref_genome <directory to reference genome> \
--in_fastq <directory to input files> \
--outfile <directory to output files> \
--sample_id <csv of sample IDs and barcode ID> \
-w <directory to save the work files>


Mandatory arguments:
      --in_fastq                  Path to input fastq dirctory 
      --outfile                   Path to output directory
      --ref_genome                reference sequence
      --sample_id                 a csv file containing barcode_Ids and sample Ids

Other arguments:
      --barcods                   barcods used during sequencing. The default 
                                  barcoding kits are "EXP-NBD104 EXP-NBD114"
      --cpu                       cpus to used during the analyis. default is 8
      --lowerlength               set the lower length for input reads filter (default: 1000)
      --upperlength               set the upper length for input reads filter (default: 20000)
      --cpu                       CPUs to be used during the analysis. The default is 8
      --memory                    Memory allocated for each process. The default is 30 GB
      --lowerlength               Set the lower length for input reads filter (default: 1000)
      --upperlength               Set the upper length for input reads filter (default: 10000)
      --pipeline                  Specify whether you want to do genome assembly or generate haplotype. The default is assembly
      --genomesize                Only required if you are running genome assembly (default: 5k)
      --basecallers               you should specify the basecalling tool you want to use with ddorado the default if basecaller and the alternative is duplex
      --model                     Please specify the spped to run basecalling, the default is sup, the alternatives are fast, hac, for more information vist dorado github
      --basecalling               Please specify whether you would like to carry out basecalling ON or OFF the default is OFF. If the basecalling is ON make sure you provide the raw fast5 or POD5 file
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
    path("basecalled/fastq_pass"), emit: in_fastq
 
    script:
    """
    #!/usr/bin/env bash
    
    # Create the output directory
    mkdir -p basecalled
    
    # Loop through input files and run dorado
    for file in \$(ls $raw_file); do
        sample_id=\$(basename "\$file")
        
        # Run the dorado command
        dorado "$basecallers" --emit-fastq "$model" "$raw_file" > "basecalled/\$sample_id"
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
* Run Barcoding
*/
process runBarcoding {
    tag "barcodes"
    publishDir params.outfile, mode: "copy", overwrite: false
            
    input:
    path (choped) 
    val (barcods)
    val (cpu)
    
    output:
    path(demultiplexed_dir), emit: demultiplexed
    path("demultiplexed_dir/*/*.fastq")

    script: 
    """
    mkdir demultiplexed_dir
    guppy_barcoder -i "${choped}" -s demultiplexed_dir \
    --barcode_kits "${barcods}" -t "${cpu}" --fastq_out --min_score_rear_override 75 --min_score 75 --trim_barcodes
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
    path(demultiplexed)
    file(ref_genome)
    val(lowerlength)
    val(upperlength)

    output:
    path("mapped_reads"), emit: mapped

    script:
    """
    #!/usr/bin/env python

    import os
    import shutil
    import subprocess

    demultiplexed_dir = "${demultiplexed}"
    ref_genome = "${ref_genome}"  # Use the provided reference genome value
    lowerlength = "${lowerlength}"
    upperlength = "${upperlength}"

    output_dir = os.path.join(os.path.dirname(demultiplexed_dir), "mapped_reads")
    os.makedirs(output_dir, exist_ok=True)

    # Copy the reference genome to the demultiplexed folder
    ref_genome_dest = os.path.join(demultiplexed_dir, os.path.basename(ref_genome))
    shutil.copyfile(ref_genome, ref_genome_dest)

    for dir in os.listdir(demultiplexed_dir):
        if not dir.startswith("bar"):
            continue

        barcode_id = dir
        dir_path = os.path.join(demultiplexed_dir, dir)

        if dir == 'barcoding_summary.txt':  # Skip the barcoding_summary.txt file
            continue

        output_subdir = os.path.join(output_dir, dir)
        os.makedirs(output_subdir, exist_ok=True)

        # Run the mapping command
        command = f"cat {dir_path}/*.fastq | \
                    minimap2 -a {ref_genome_dest} /dev/stdin | \
                    samtools view -b | \
                    samtools sort > {output_subdir}/{barcode_id}_align_sorted.bam"
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
* Run error_correction
*/

process runErrcorrect {
    tag {"error_correction"}
    label 'small_mem'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(mapped)

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
    subdirs = [d for d in os.listdir(mapped_dir) if os.path.isdir(os.path.join(mapped_dir, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(mapped_dir, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Run the error correction commands
        command1 = f"${baseDir}/packages/RATTLE/rattle cluster -i {subdir_path}/{subdir}_mapped_reads.fastq -p 0.2 -t 8 --iso --repr-percentile 0.3 -o {output_subdir}"
        subprocess.run(command1, shell=True, check=True)

        command2 = f"${baseDir}/packages/RATTLE/rattle correct -i {subdir_path}/{subdir}_mapped_reads.fastq -c {output_subdir}/clusters.out -t 8 -o {output_subdir}"
        subprocess.run(command2, shell=True, check=True)

        # Move the corrected file
        corrected_file = os.path.join(output_subdir, f"{subdir}_corrected.fastq")
        os.makedirs(os.path.dirname(corrected_file), exist_ok=True)
        command3 = f"mv {output_subdir}/corrected.fq {corrected_file}"
        subprocess.run(command3, shell=True, check=True)
    """
}

/*
* Run genome_assembly
*/

process runAssembly {
    tag {"genome_assembly"}
    label 'bonobo_img'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(corrected)
    val(genomesize)
    
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

        # Run the assembly command
        command1 = f"flye --meta --genome-size {genomesize} --nano-raw {subdir_path}/{subdir}_corrected.fastq --no-alt-contigs --out-dir {output_subdir}"
        subprocess.run(command1, shell=True, check=True)
        
        # Move the assembly file
        assembled_file = os.path.join(output_subdir, f"{subdir}_contigs.fasta")
        os.makedirs(os.path.dirname(assembled_file), exist_ok=True)
        command2 = f"mv {output_subdir}/assembly.fasta {assembled_file}"
        subprocess.run(command2, shell=True, check=True)

    """
}

/*
* Run haplotype
*/

process runHaplotype {
    tag {"genome_haplotype"}
    label 'bonobo_img'
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(corrected)
    
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

        # Run the seqtk command
        command1 = f"seqtk seq -A {subdir_path}/{subdir}_corrected.fastq > {output_subdir}/{subdir}_corrected.fasta"
        subprocess.run(command1, shell=True, check=True)

        # Run the strainline command
        command2 = f"/app/Strainline/src/strainline.sh -i {output_subdir}/{subdir}_corrected.fasta -o {output_subdir} -p ont --maxLD 0.01 --rmMisassembly True --threads ${params.cpu}"
        subprocess.run(command2, shell=True, check=True)
        
        # Move the haplotype file
        haplotype_file = os.path.join(output_subdir, f"{subdir}_contigs.fasta")
        os.makedirs(os.path.dirname(haplotype_file), exist_ok=True)
        command3 = f"mv {output_subdir}/haplotypes.final.fa {haplotype_file}"
        subprocess.run(command3, shell=True, check=True)
    """
}


/*
* Run polishing_medaka
*/

process runPolish_medaka {
    tag {"polish_medaka"}
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(draftgenome)
    path(mapped)
    
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

        # Print the file paths for debugging
        print("Assembled contigs file:", os.path.join(subdir_path, f"{subdir}_contigs.fasta"))
        print("Mapped reads file:", mapped_file)

        # Check the file sizes for debugging
        print("Assembled contigs file size:", os.path.getsize(os.path.join(subdir_path, f"{subdir}_contigs.fasta")))
        print("Mapped reads file size:", os.path.getsize(mapped_file))

        # Run the medaka command
        command1 = f"medaka_consensus -i {mapped_file} \
                         -d {subdir_path}/{subdir}_contigs.fasta \
                         -o {output_subdir} \
                         -t 8"
        subprocess.run(command1, shell=True, check=True)
        
        # Move the medaka file
        polish_medaka_file = os.path.join(output_subdir, f"{subdir}_polish_medaka.fasta")
        os.makedirs(os.path.dirname(polish_medaka_file), exist_ok=True)
        command2 = f"mv {output_subdir}/consensus.fasta {polish_medaka_file}"
        subprocess.run(command2, shell=True, check=True)
    """
}

/*
* Run polishing_pilon
*/

process runPolish_pilon {
    tag {"polish_pilon"}
    publishDir "${params.outfile}", mode: "copy", overwrite: false

    input:
    path(medaka)
    path(mapped)
    
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

    # Get the subdirectories in the assembled directory
    subdirs = [d for d in os.listdir(polished) if os.path.isdir(os.path.join(polished, d))]

    for subdir in subdirs:
        if not subdir.startswith("bar"):
            continue

        subdir_path = os.path.join(polished, subdir)
        output_subdir = os.path.join(output_dir, subdir)
        os.makedirs(output_subdir, exist_ok=True)

        # Get the corresponding mapped file
        mapped_file = os.path.join(mapped_dir, subdir, f"{subdir}_mapped_reads.fastq")

        # Print the file paths for debugging
        print("Assembled contigs file:", os.path.join(subdir_path, f"{subdir}_polish_medaka.fasta"))
        print("Mapped reads file:", mapped_file)

        # Check the file sizes for debugging
        print("Assembled contigs file size:", os.path.getsize(os.path.join(subdir_path, f"{subdir}_polish_medaka.fasta")))
        print("Mapped reads file size:", os.path.getsize(mapped_file))

        # Run the pilon command
        command1 = f"bwa index {subdir_path}/{subdir}_polish_medaka.fasta && \
                     bwa mem -t7 {subdir_path}/{subdir}_polish_medaka.fasta {mapped_file} | \
                     samtools sort -@ 7 | \
                     samtools markdup /dev/stdin {output_subdir}/{subdir}_cons_ali_filter.bam && \
                     samtools view -b -@ 7 -q 30 {output_subdir}/{subdir}_cons_ali_filter.bam \
                                          -o {output_subdir}/{subdir}_cons_ali_dup_filter.bam && \
                     samtools index -@ 7 {output_subdir}/{subdir}_cons_ali_dup_filter.bam && \
                     pilon -Xmx4096m -XX:-UseGCOverheadLimit \
                           --genome {subdir_path}/{subdir}_polish_medaka.fasta \
                           --unpaired {output_subdir}/{subdir}_cons_ali_dup_filter.bam \
                           --output {output_subdir}/{subdir}_polishing \
                           --vcf"
        subprocess.run(command1, shell=True, check=True)
    """
}


/*
* Run polishing_prooveframe
*/


process runProovframe {
    label 'bonobo_img'
    tag {"proovframe"}
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

        # Run the proovframe commands
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

    # Read the CSV file and create the mapping
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


workflow {
    if (params.basecalling == 'OFF') {
        runPowerchop(params.in_fastq, params.cpu)
    }
    else if (params.basecalling == 'ON') {
        runDorado(params.raw_file, params.cpu, params.basecallers, params.model)
        runPowerchop(runDorado.out.in_fastq, params.cpu)
        }

   runBarcoding(runPowerchop.out.choped, params.barcods, params.cpu)
   runMapping(runBarcoding.out.demultiplexed, params.ref_genome, params.lowerlength, params.upperlength)
   runErrcorrect(runMapping.out.mapped)

    if (params.pipeline == 'assembly') {
       runAssembly(runErrcorrect.out.corrected, params.genomesize)
       runPolish_medaka(runAssembly.out.draftgenome, runMapping.out.mapped)
   }

    else if (params.pipeline == 'haplotype') {
       runHaplotype(runErrcorrect.out.corrected)
       runPolish_medaka(runHaplotype.out.draftgenome, runMapping.out.mapped)
   }

   runPolish_pilon(runPolish_medaka.out.medaka, runMapping.out.mapped)
   runProovframe(runPolish_pilon.out.pilon)
   runSeqrenaming(runProovframe.out.final_seq, params.sample_id)
}