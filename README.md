# NanoRepeat

NanoRepeat is a computational tool for tandem repeat quantification from Oxford Nanopore long-read sequencing data

## Table of Contents

- [Workflow](#Workflow)
- [Prerequisites](#Prerequisites)
- [Installation](#Installation)
- [Usage](#Usage)


## <a name="Workflow"></a>1 Workflow

NanoRepeat can quantify a single repeat, or perform a joint quantification of two adjacent repeats (such as the `CAG` and `CCG` repeats in the HTT gene).

### 1.1 Quantification of a single repeat (most common)
NanoRepeat generates a series of sequences where the repeat sizes are from 1 to N (a user specified value) with 10 kb left and right flanking sequences. The reads were aligned to this series of sequences using [minimap2](https://github.com/lh3/minimap2) with the parameter for the specified platform. The repeat size of the sequence with the highest alignment score was the estimate of the repeat size of the read. 

<p align="center"><img src="images/repeat_size_estimation.jpg" alt="repeat_size_estimation" width="100%"></p>

After the repeat size of each read is determined, NanoRepeat uses the Gaussian mixture model (GMM) to assign reads to alleles. First, outlier reads with repeat sizes that are outside three standard deviations from the mean are removed. Next, we assume that the repeat size is a mixture of N Gaussian models, where N is 1 or 2 for a diploid genome. We use the Bayesian information criterion (BIC) to select the best N. After the best N is selected, the label of each read is predicted using the trained Gaussian mixture model. In high-confidence mode, a read is discarded if it is not confidently assigned to a Gaussian model (p < 0.95) or if it is outside three standard deviations from the mean of that model.

<p align="center"><img src="images/classify_reads.jpg" alt="classify_reads" width="100%"></p>

### 1.2 Joint quantification two adjacent repeats

In some regions, there are two adjacent repeats with very similar sequences. For example, in the HTT gene (cause Huntington's disease), the repeat is `(CAG)m-CAA-CAG-CCG-CCA-(CCG)n` and the CAG and CCG repeats are of variable sizes. In this case, NanoRepeat can perform a joint quantification of the two repeats. 

![image](https://user-images.githubusercontent.com/9782948/142343820-94211ca0-8f25-4a1f-8991-332b7659c7b4.png)

The joint quantification process has two steps: fast estimation and refining. In the fast estimation step, NanoRepeat performs a quick analysis of each repeat and estimates the lower and upper bound of the two repeat sizes. Let L1, L2 denote the lower bounds of the CAG and CCG repeats, and U1, U2 denote the upper bounds of the two repeats, respectively. In the refining step, NanoRepeat generates a batch of amplicon sequences with m (L1 ≤ m ≤ U1) CAG repeat units and n (L2 ≤ n ≤ U2) CCG repeat units. Each read is aligned to this batch of amplicon sequences with minimap2. The m and n that maximize the alignment score were the estimated CAG and CCG repeat sizes of the read.

After the repeat number of each read was determined, NanoRepeat classifies the reads to alleles using a 2D GMM model. Outlier reads are removed as shown in the following figure. 

![image](https://user-images.githubusercontent.com/9782948/142344767-f41787f1-29a4-4f53-821e-e9b2d01ab797.png)



## <a name="Prerequisites"></a>2 Prerequisites

1. Operating system: Linux or MacOS
2. Python3 (Python 2 is NOT supported)
3. Python packages: `numpy`, `sklearn`, `matplotlib`. You can use `pip` to install the packages:
    ```
    pip install --user numpy sklearn matplotlib
    ```
4. Minimap2 (version >= 2.8). NanoRepeat calls `minimap2` to do sequence alignment. If you don't have `minimap2` in your system, you can install it following the instructions [here](https://github.com/lh3/minimap2#install). If you are using Linux, you can acquire precompiled binaries using the following commands:

    ```
    wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2
    tar -jxvf minimap2-2.24_x64-linux.tar.bz2
    ```

5. Samtools (version >= 1.3). NanoRepeat calls `samtools` to process SAM and BAM files. If you don't have `samtools` in your system, you can download the source code and compile it.

    ```
    wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
    tar -jxvf samtools-1.7.tar.bz2
    cd samtools-1.7
    make
    ```
   Please note that some libraries (e.g. HTSlib) are needed to compile `samtools`. If you are unable to compile `samtools` from the source code, you can install `samtools` using Anaconda: `conda install -c bioconda samtools`. 
    
## <a name="Installation"></a>3 Installation

You can clone the repository of NanoRepeat using the following command.
```
git clone https://github.com/WGLab/NanoRepeat.git
```

The scripts in the `./NanoRepeat` can run directly without additional compilation or installation.

## <a name="Usage"></a>4 Usage

```
usage: nanoRepeat.py [-h] -i input_file -t input_type -r ref.fasta -b repeat_regions.bed -o prefix/of/output/files [-c INT] [--samtools path/to/samtools] [--minimap2 path/to/minimap2] [--ploidy INT] [--anchor_len INT]

NanoRepeat: short tandem repeat (STR) quantification from Nanopore long-read sequencing

optional arguments:
  -h, --help            show this help message and exit
  -i input_file, --input input_file
                        (required) path to input file (supported format: sorted_bam, fastq or fasta)
  -t input_type, --type input_type
                        (required) input file type (valid values: bam, fastq or fasta)
  -r ref.fasta, --ref_fasta ref.fasta
                        (required) path to reference genome sequence in FASTA format
  -b repeat_regions.bed, --repeat_region_bed repeat_regions.bed
                        (required) path to repeat region file (in bed format)
  -o prefix/of/output/files, --out_prefix prefix/of/output/files
                        (required) prefix of output files
  -c INT, --num_cpu INT
                        (optional) number of CPU cores (default: 1)
  --samtools path/to/samtools
                        (optional) path to samtools (default: using environment default)
  --minimap2 path/to/minimap2
                        (optional) path to minimap2 (default: using environment default)
  --ploidy INT          (optional) ploidy of the sample (default: 2)
  --anchor_len INT      (optional) length of up/downstream sequence to help identify the repeat region (default: 256 bp, increase this value if the 1000 bp up/downstream sequences are also repeat)

Examples: 
	1) python nanoRepeat.py -i input.bam   -t bam   -r hg38.fasta -b hg38.repeats.bed -c 4 -o prefix/of/output/files
	2) python nanoRepeat.py -i input.fastq -t fastq -r hg38.fasta -b hg38.repeats.bed -c 4 -o prefix/of/output/files
	3) python nanoRepeat.py -i input.fasta -t fasta -r hg38.fasta -b hg38.repeats.bed -c 4 -o prefix/of/output/files
```
