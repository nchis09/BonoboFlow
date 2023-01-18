#  BonoboFlow
[![](https://img.shields.io/badge/nextflow-22.10.4-yellowgreen)](https://www.nextflow.io)
[![](https://img.shields.io/badge/uses-docker-orange)](https://docs.docker.com/get-docker)
[![](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[![Twitter Follow](https://img.shields.io/twitter/follow/ndekezi09.svg?style=social)](https://twitter.com/ndekezi09) 


##                      B O N O B O F L O W - P I P E L I N E


     BonoboFlow PIPELINE for viral genome assemblypipeline from MinION 

                                  sequenced reads

                ------------------------------------------------------

                      ------------------------------------------
                
                                --------------------
                
                                        -----
                            
                                          -

                      =============================================
                         BonoboFlow  ~  version 1.0

                                  Athours: Ndekezi Christian
                      ==============================================


## Introduction

BonoboFlow is a nextflow pipeline for reproducible analysis of 



## Installation

BonoboFlow requires **nextflow** (version 22.10.4), docker, Diamond (version  2.9.15 or higher), minimap2, samtools, bamtools, filtlong)

conda create -n nextflow -c bioconda -c conda-forge openjdk=11.0.8 nextflow
for mac uses we recommend to install diamond using homebrew "brew install diamond"
pip install isONcorrect


Note: make sure docker has at least 60GB available disk space, 6CPUS, 16GB memory


```bash
git clone 
nextflow run bonoboflow.nf --help
```

## Usage

The pipeline takes in the fast5 file from Nanopore sequencing technology 

```bash

nextflow run BonoboFlow_2.nf -resume --kit <sequencing kit>\
           --flowcell <flow cell used during sequencing> \
           --cpu allocated_cpu --in_fast5 <directory to input files> \
           --outfile <directory to output files> \
           --ref_genome <directory to reference genome>
    

    Mandatory arguments:
      --kit                       Sequencit kit 
      --flowcell                  Flowcell used during sequencing
      --in_fast5                  Path to input fast5 dirctory 
      --outfile                   Path to output directory
      --ref_genome                reference sequence

    Other arguments:
      --barcods                   barcods used during sequencing. The default barcoding kits are "EXP-NBD104 EXP-NBD114"
      --cpu                       cpus to used during the analyis. default is 8
      --lowerlength               set the lower length for input reads filter (default: 150)

```

### Profiles

BonoboFlow can be run under different computing environments, simply choose an appropriate profile via the `-profile` argument. Could take only `-profile docker`


#### Pipeline 


**Output fiiles**

* `consensus*.fasta`: `FASTA` files of consensus sequences - one per sample


**Pipeline information output**

* `QuasiFlow_DAG.html`: Graphical representation of the pipeline's processes/operators and channels between them.



## Parameters

Mandatory parameters

* `--kit `: Sequencit kit 
* `--flowcell`:              Flowcell used during sequencing
* `--in_fast5`:                Path to input fast5 dirctory 
* `--ref_genome`:                reference sequence


Optional parameters

* `--barcods`:        barcods used during sequencing. The default barcoding kits are "EXP-NBD104 EXP-NBD114"
* `--cpu`:                 cpus to used during the analyis. default is 8
* `--lowerlength`:               set the lower length for input reads filter (default: 150)
* `--upperlength`:             set the upper length for input reads filter (default: 100,000)


### Output parameters

Mandatory parameters

* `--outfile`:          Path to output directory

## Dependancies.

Below is the list of tools that are used in the BonoboFlow pipeline. These tools are readliy available and may be installed using `conda` via `bioconda` channel.

+ [minimap2](https://github.com/lh3/minimap2)
+ [samtools](https://github.com/samtools/samtools)
+ [bamtools](https://github.com/pezmaster31/bamtools)
+ [diamond](https://github.com/bbuchfink/diamond)




## Troubleshooting

Kindly report any issues at 

## License

BonoboFlow is licensed under GNU GPL v3.

## Citation

**This work is currently under peer review. A formal citation will be availed in due course.**

# BonoboFlow
