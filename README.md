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
conda create -n nextflow -c bioconda -c conda-forge openjdk=11.0.8 nextflow python cmake spoa && \
conda activate nextflow && \
cd BonoboFlow/packages/RATTLE && \
./build.sh && \
cd ../../
```

Note: make sure docker has at least 60GB available disk space, 6CPUS, 16GB memory


```bash

nextflow run BonoboFlow.nf --help
```

## Usage

The pipeline takes in the fast5 file from Nanopore sequencing technology 

```bash
conda activate nextflow

nextflow run BonoboFlow.nf -resume --ref_genome <directory to reference genome> --in_fastq <directory to input files> --outfile <directory to output files> --sample_id <csv of sample IDs and barcode ID> -w <directory to save the work files>

    

Mandatory arguments:
      --in_fastq                  Path to input fastq dirctory. Note: If you specify this you dont have to specify the  --raw_file
      --raw_file                  Path to raw POD5 or FAST5 files. Note: If you specify this, make sure you change the --basecalling flag to ON
      --outfile                   Path to output directory
      --ref_genome                reference sequence
      --sample_id                 a csv file containing barcode_Ids and sample Ids. An example csv file is provided in the BonoboFlow directory

Other arguments:
      --barcods                   barcods used during sequencing. The default barcoding kits are "EXP-NBD104 EXP-NBD114"
      --cpu                       cpus to used during the analyis. default is 8
      --lowerlength               set the lower length for input reads filter (default: 1000)
      --upperlength               set the upper length for input reads filter (default: 20000)
      --memory                    Memory allocated for each process. The default is 30 GB
      --pipeline                  Specify whether you want to do genome assembly or haplotype reconstruction. The default is haplotype
      --genomesize                Only required if you are running genome assembly (default: 5k)
      --basecalling               Please specify whether you would like to carry out basecalling ON or OFF the default is OFF. If the basecalling is ON make sure you provide the raw fast5 or POD5 file
      --basecallers               you should specify the basecalling tool you want to use with ddorado the default if basecaller and the alternative is duplex
      --model                     Please specify the spped to run basecalling, the default is sup, the alternatives are fast, hac, for more information vist dorado github



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
