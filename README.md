# BonoboFlow
[![](https://img.shields.io/badge/nextflow-24.04.2-yellowgreen)](https://www.nextflow.io)
[![](https://img.shields.io/badge/uses-docker-orange)](https://docs.docker.com/get-docker)
[![](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![Twitter Follow](https://img.shields.io/twitter/follow/ndekezi09.svg?style=social)](https://twitter.com/ndekezi09) 

## Introduction

BonoboFlow is a Nextflow pipeline for reproducible and precise execution of viral genome assembly and haplotypes construction from MinION sequencing reads. The pipeline supports both raw FAST5/POD5 files and FASTQ inputs, with capabilities for basecalling, demultiplexing, error correction, assembly, and haplotype reconstruction.

## Features

- Basecalling from FAST5/POD5 files using Dorado
- Demultiplexing with support for multiple barcoding kits
- Read quality filtering and length selection
- Error correction using either Rattle or Vechat
- Genome assembly or haplotype reconstruction
- Automatic container detection and configuration (Docker/Singularity)
- Flexible resource allocation for different computational requirements
- DAG visualization for workflow monitoring

## Installation

### Prerequisites
- Nextflow (version 24.04.2)
- Docker or Singularity
- Conda (recommended for environment management)

### Quick Start
```bash
# Clone the repository and create conda environment
git clone https://github.com/nchis09/BonoboFlow.git && \
conda create -n bonoboflow -c bioconda -c conda-forge openjdk=11 nextflow=24.04.2 python && \
conda activate bonoboflow

# Test installation
conda activate bonoboflow
nextflow run BonoboFlow.nf --help
```

## Container Support

BonoboFlow automatically detects and uses either Docker or Singularity:
- Docker: Runs with user permissions using `-u $(id -u):$(id -g)`
- Singularity: Runs Docker images
- Containers are pulled automatically for each process
- Default container images:
  - Main pipeline: `nchis09/bonobo_image:v2`
  - Basecalling/Barcoding: `nanoporetech/dorado:latest`
  - Assembly corrections: `nanoporetech/medaka:latest`
  - Prooveframe: `nanoporetech/prooveframe:latest`

## Resource Configuration

### Default Settings
- CPU cores: 8
- Memory: 16 GB

### Process-Specific Resources
- Main pipeline processes: Up to 300 GB memory with automatic retry on memory errors
- Basecalling and demultiplexing: Configurable based on input parameters

## Usage

### Basic Command Structure
```bash
nextflow run BonoboFlow.nf -resume \
    --in_fastq <input_directory> \
    --outfile <output_directory> \
    --ref_genome <reference_genome> \
    --sample_id <sample_csv_file> \
    -w <work_directory>
```

### Input Parameters

#### Mandatory Arguments
| Parameter | Description |
|-----------|-------------|
| `--in_fastq` | Path to input FASTQ directory (mutually exclusive with --raw_file) |
| `--raw_file` | Path to raw POD5/FAST5 files (requires --basecalling ON) |
| `--outfile` | Output directory path |
| `--ref_genome` | Reference genome sequence |
| `--sample_id` | CSV file with barcode and sample IDs |

#### Optional Arguments
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cpu` | 8 | Number of CPUs to use |
| `--memory` | "16 GB" | Memory allocation per process |
| `--pipeline` | "assembly" | Pipeline mode: "assembly" or "haplotype" |
| `--phred` | 6 | Minimum sequence quality score |
| `--lowerlength` | 1000 | Minimum read length |
| `--upperlength` | 20000 | Maximum read length |
| `--genomesize` | "5k" | Expected genome size (assembly mode only) |

#### Performance Tuning
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cpu` | 50 | Number of CPUs to use |
| `--memory` | "50 GB" | Memory allocation per process |
| `--split_size` | 10000 | Split size for parallel processing |

#### Quality Control
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--phred` | 6 | Minimum sequence quality score |
| `--lowerlength` | 1000 | Minimum read length |
| `--upperlength` | 200000 | Maximum read length |
| `--min_mq` | 20 | Minimum mapping quality |

#### Basecalling Configuration
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--basecalling` | "OFF" | Enable/disable basecalling |
| `--basecallers` | "basecaller" | Tool choice: "basecaller" or "duplex" |
| `--model` | "sup" | Model type: "sup", "fast", or "hac" |

#### Barcoding Options
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--barcods` | "SQK-NBD114-96" | Barcoding kits used |
| `--min_score_rear_barcode` | 75 | Minimum rear barcode quality |
| `--min_score_front_barcode` | 75 | Minimum front barcode quality |

#### Error Correction Settings
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--error_correction_tool` | "vechat" | Choose between "vechat" or "rattle" |
| `--repr_percentile` | 0.3 | Cluster representative percentile |
| `--score_threshold` | 0.2 | Cluster similarity threshold |
| `--kmer_size` | 12 | K-mer size for clustering |

#### Haplotype Construction
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--maxLD_floats` | 0.001 | Maximum local divergence |
| `--maxGD_floats` | 0.01 | Maximum global divergence |
| `--minAbun_floats` | 0.01 | Minimum abundance threshold |
| `--rmMisassembly_bool` | True | Remove misassemblies |
| `--correctErr` | False | Error correction for input reads |
| `--topks` | 200 | Number of seed reads for haplotype construction (higher value = more seed reads) |
| `--minovlplens` | 500 | Minimum read overlap length (lower value = more permissive overlaps) |
| `--maxohs` | 30 | Maximum overhang length for overlaps (higher value = more permissive overlaps) |

> **Note**: Error correction is now disabled by default in the latest version to improve processing speed. Enable it with `--correctErr True` if needed for your specific use case.

> **Note**: These parameters are optimized for maximum haplotype detection. They are more permissive than the default Strainline settings to capture more potential haplotypes, including rare variants. Validate results carefully as these settings might increase false positives.

## Output Files

- `consensus*.fasta`: Final consensus sequences per sample
- `final_reports/`: Directory containing:
  - `BonoboFlow_timeline.html`: Process execution timeline
  - `BonoboFlow_report.html`: Detailed execution report
  - `BonoboFlow_DAG.html`: Pipeline workflow visualization
  - `consensus*.fasta`: Final consensus sequences per sample


## Troubleshooting

For issues and support:
1. Check the execution reports in `final_reports/`
2. Submit issues at https://github.com/nchis09/BonoboFlow/issues
3. Include relevant error messages and configuration

## License

This project is licensed under the GNU General Public License v3.0.

## Citation

**Christian Ndekezi, Drake Byamukama, Frank Kato, Denis Omara, et al, BonoboFlow: Viral Genome Assembly and Haplotype Reconstruction from Nanopore Reads, Bioinformatics Advances, 2025;, vbaf115, https://doi.org/10.1093/bioadv/vbaf115**