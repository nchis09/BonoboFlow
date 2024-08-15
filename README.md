#  BonoboFlow
[![](https://img.shields.io/badge/nextflow-22.10.4-yellowgreen)](https://www.nextflow.io)
[![](https://img.shields.io/badge/uses-docker-orange)](https://docs.docker.com/get-docker)
[![](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[![Twitter Follow](https://img.shields.io/twitter/follow/ndekezi09.svg?style=social)](https://twitter.com/ndekezi09) 


## B O N O B O F L O W - P I P E L I N E


The BonoboFlow pipeline is a dedicated tool developed for the precise execution of viral 
        genome assembly and haplotypes construction from MinION sequencing reads

                ------------------------------------------------------

                      ------------------------------------------
                
                                --------------------
                
                                        -----
                            
                                          -

                      =============================================
                              BonoboFlow  ~  version 1.0
                      ==============================================


## Introduction

BonoboFlow is a nextflow pipeline for reproducible and precise execution of viral genome assembly and haplotypes construction 


## Installation

BonoboFlow requires:
 **nextflow** (version 22.10.4)
 docker
 singularity


We recommend to install the packages as follow

```bash
git clone https://github.com/nchis09/BonoboFlow.git && \
conda create -n bonoboflow -c bioconda -c conda-forge openjdk=11.0.8 nextflow=24.04.2  python cmake spoa && \
conda activate bonoboflow && \
cd BonoboFlow/packages/RATTLE && \
./build.sh && \
cd ../../
```

```bash
conda activate bonoboflow
nextflow run BonoboFlow.nf --help
```

## Usage

The pipeline takes in the fast5 file from Nanopore sequencing technology 

```bash
conda activate bonoboflow

nextflow run BonoboFlow.nf -resume --ref_genome <directory to reference genome> --in_fastq <directory to input files> --outfile <directory to output files> --sample_id <csv of sample IDs and barcode ID> -w <directory to save the work files>

    

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

Error_correction arguments:
      --repr_percentile           cluster representative percentile (default: 0.15)
      --score_threshold           minimum score for two reads to be in the same gene cluster (default: 0.2)
      --kmer_size                 k-mer size for isoform clustering (default: 11, maximum: 16)

Haplotype arguments:
      --maxLD_floats              Maximum local divergence allowed for merging haplotypes. (default: 0.01)
      --maxGD_floats              Maximum global divergence allowed for merging haplotypes. (default: 0.01)
      --rmMisassembly_bool        Break contigs at potential misassembled positions (default: False)
      --correctErr_bool           Perform error correction for input reads (default: False)
      --minAbun_floats            Minimum abundance for filtering haplotypes (default: 0.02)
      --topks                     k seed reads size for haplotype construction (default: 100)
      --minovlplens               Minimum read overlap length. (default: 1000)
      --minseedlens               Minimum seed read length. (default: 2000)
      --maxohs                    Maximum overhang length allowed for overlaps. (default: 20)


```

### Profiles

BonoboFlow can be run under different computing environments, simply choose an appropriate profile via the `-profile` argument. Could take only `-profile docker`


#### Pipeline 


**Output fiiles**

* `consensus*.fasta`: `FASTA` files of consensus sequences - one per sample


**Pipeline information output**

* `Bonoboflow_DAG.html`: Graphical representation of the pipeline's processes/operators and channels between them.



## input Parameters

Mandatory parameters

* `--in_fastq `:            Path to input fastq dirctory
* `--sample_id`:           a csv file containing barcode_Ids and sample Ids
* `--ref_genome`:           reference sequence


Optional parameters

* `--barcods`:        barcods used during sequencing. The default barcoding kits are "EXP-NBD104 EXP-NBD114"
* `--cpu`:                 cpus to used during the analyis. default is 8
* `--lowerlength`:               set the lower length for input reads filter (default: 150)
* `--upperlength`:             set the upper length for input reads filter (default: 100,000)


### Output parameters

Mandatory parameters

* `--outfile`:          Path to output directory


## Troubleshooting

Kindly report any issues at https://github.com/nchis09/BonoboFlow/issues

## License

BonoboFlow is licensed under GNU GPL v3.

## Citation

**This work is currently under peer review. A formal citation will be availed in due course.**

# BonoboFlow
