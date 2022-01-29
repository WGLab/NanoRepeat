# NanoRepeat

NanoRepeat is a computational tool for tandem repeat detection from Nanopore long-read sequencing data. 

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
    wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
    tar -jxvf minimap2-2.17_x64-linux.tar.bz2
    ./minimap2-2.17_x64-linux/minimap2
    ```
## <a name="Installation"></a>3 Installation

You can clone the repository of NanoRepeat using the following command.
```
git clone https://github.com/WGLab/NanoRepeat.git
```

The scripts in the `./NanoRepeat` can run directly without additional compilation or installation.

## <a name="Usage"></a>4 Usage

```
usage: NanoRepeat.py [-h] --in_fq PATH --platform STRING --ref_fasta PATH
                    --repeat_region chr:start-end --repeat_unit STRING
                    --max_repeat_size INT --out_dir PATH [--num_threads INT]
                    [--minimap2 PATH] [--fixed_cutoff_value INT]
                    [--ploidy INT] [--anchor_len INT] [--version]

A tools for short tandem repeat detection from long-read amplicon sequencing
data

optional arguments:
  -h, --help            show this help message and exit
  --in_fq PATH          path to input fastq file
  --platform STRING     three valid values: ont, pacbio, consensus
  --ref_fasta PATH      path to reference genome sequence in FASTA format
  --repeat_region chr:start-end
                        repeat region in the reference genome (e.g.
                        chr4:3074876-3074939, coordinates are 1-based)
  --repeat_unit STRING  sequence of the repeat unit (e.g. CAG)
  --max_repeat_size INT
                        maximum possible number of the repeat unit
  --out_dir PATH        path to the output directory
  --num_threads INT     number of threads used by minimap2 (default: 1)
  --minimap2 PATH       path to minimap2 (default: using environment default)
  --fixed_cutoff_value INT
                        split alleles using this fixed_cutoff_value (if set,
                        ploidy will be set to 2). reads with repeat size >=
                        fixed_cutoff_value will be assigned to the second
                        allele
  --ploidy INT          ploidy of the sample (default: 2)
  --anchor_len INT      length of up/downstream sequence to help identify the
                        repeat region (default: 1000 bp, increase this value
                        if the 1000 bp up/downstream sequences are also
                        repeat)
  --version             show program's version number and exit
```
