#  BonoboFlow
[![](https://img.shields.io/badge/nextflow-22.10.4-yellowgreen)](https://www.nextflow.io)
[![](https://img.shields.io/badge/uses-docker-orange)](https://docs.docker.com/get-docker)
[![](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[![Twitter Follow](https://img.shields.io/twitter/follow/ndekezi09.svg?style=social)](https://twitter.com/ndekezi09) 


## B O N O B O F L O W - P I P E L I N E


        BonoboFlow PIPELINE for viral genome assemblypipeline from MinION 

                                  sequenced reads

                ------------------------------------------------------

                      ------------------------------------------
                
                                --------------------
                
                                        -----
                            
                                          -

                      =============================================
                         BonoboFlow  ~  version 1.0
                      ==============================================


## Introduction

BonoboFlow is a nextflow pipeline for reproducible analysis of 



## Installation

BonoboFlow requires:
 **nextflow** (version 22.10.4)
 docker
 singularity




We recommend to install the packages as follow

```bash
conda create -n nextflow -c bioconda -c conda-forge openjdk=11.0.8 nextflow python
conda activate nextflow
```

Note: make sure docker has at least 60GB available disk space, 6CPUS, 16GB memory


```bash
git clone https://github.com/nchis09/BonoboFlow.git
nextflow run bonoboflow.nf --help
```

## Usage

The pipeline takes in the fast5 file from Nanopore sequencing technology 

```bash
conda activate nextflow

nextflow run BonoboFlow.nf -resume --in_fastq <in put directory> --outfile <output directory> --ref_genome <reference genone.fasta> --lowerlength 1000 --upperlength 5000 -w <work directory> --sample_id <sample_id.csv> --kit <sequencing kit>  --flowcell <flow cell used during sequencing> 

    

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
