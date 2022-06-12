# NanoRepeat: quantification of Short Tandem Repeats (STRs) from Oxford Nanopore long-read sequencing data

## Table of Contents

- [Installation](#installation)
- [Repeat quantification](#repeat_quantification)
  - [Joint quantification of two adjacent repeats](#joint_quantification)


## <a name="installation"></a> Installation

NanoRepeat requires the following software tools:

1. [Python](https://www.python.org/downloads/) (version >= 3.7)
2. [Minimap2](https://github.com/lh3/minimap2) (version >= 2.8)
3. [Samtools](https://github.com/samtools/samtools.git) (version >= 1.3)

You may alreadly have `minimap2` and `samtools` if you performed analysis of Oxford Nanopore sequencing data. You can use `which minimap2` and `which samtools` to check the full path to the two executable files.

Once you installed the above tools, you can use the following commands to install NanoRepeat:
```
git clone https://github.com/WGLab/NanoRepeat.git
cd NanoRepeat
pip install -r requirements.txt
```
## <a name="repeat_quantification"></a> Repeat quantification

### <a name="joint_quantification"> Joint quantification of two adjacent STRs (such as the `CAG` and `CCG` repeats in the HTT gene)

In exon-1 of the human HTT gene, there are two adjacent STRs: `CAG` and `CCG`. The sequence structure is: (CAG)<sub>m</sub>-CAA-CAG-CCG-CCA-(CCG)<sub>n</sub>. NanoRepeat can jointly quantify the two STRs and provide phased results. In our experience, looking at both repeats help generate better quantification results. 
	
We will demonstrate the usage of NanoRepeat using an example data set. 
	
```
mkdir joint_quantification_HTT
cd joint_quantification_HTT
wget https://github.com/WGLab/NanoRepeat/releases/download/v1.0/NanoRepeat_example_data.tar.gz
tar xzf NanoRepeat_example_data.tar.gz
```
	
The input fastq file is here: `./NanoRepeat_example_data/HTT_amplicon.fastq.gz`.
	
The reference fasta file is here: `./NanoRepeat_example_data/GRCh38_chr4.0_4Mb.fasta`.

Joint quantification: 
```
python path/to/NanoRepeat/nanoRepeat-joint.py  \
    --in_fq ./NanoRepeat_example_data/HTT_amplicon.fastq.gz        \
    --ref_fasta ./NanoRepeat_example_data/GRCh38_chr4.0_4Mb.fasta  \
    --repeat1 chr4:3074876:3074933:CAG:200                         \
    --repeat2 chr4:3074946:3074966:CCG:20                          \
    --out_dir ./nanorepeat_output1                                 \
    --minimap2 path/to/minimap2                                    \ # optional if the minimap2 binary is in your environment
    --num_threads 4
```

`--repeat1` and `--repeat2` specify the two repeat regions. The format of `--repeat1`  and `--repeat2` is `chrom:start_position:end_position:repeat_unit:max_size`. The start and end positions are 0-based (the first base on the chromosome is numbered 0). The start position is self-inclusive but the end position is non-inclusive, which is the same as the [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). For example, a region of the first 100 bases of chr1 is denoted as `chr1:0:100`.  `max_size` is the max repeat length that we consider. Please set `max_size` to be a reasonal number. If `max_size` is too large (e.g. well beyond the max possible number), the speed of joint quantification might be slow. 


You will see the following files in the `./nanorepeat_output1` folder if NanoRepeat ran succesfully. 

```
HTT_amplicon.fastq.JointGMM.chr4_3074876_3074933_CAG.hist.png
HTT_amplicon.fastq.JointGMM.chr4_3074946_3074966_CCG.hist.png
HTT_amplicon.fastq.JointGMM.hist2d.png
HTT_amplicon.fastq.JointGMM.qc_passed.allele1.fastq
HTT_amplicon.fastq.JointGMM.qc_passed.allele2.fastq
HTT_amplicon.fastq.JointGMM.scatter.png
HTT_amplicon.fastq.JointGMM.summary.txt
HTT_amplicon.fastq.repeat_size.txt
```

The `*qc_passed.allele1.fastq` file is the reads assigned to the first allele and the `qc_passed.allele2.fastq` is the reads assigned to the second allele. 

The `*summary.txt` file gives the quantification of the repeat sizes. 

For example, 
```
$ cat HTT_amplicon.fastq.JointGMM.summary.txt
#input_fastq	method	num_alleles	gmm_cov_repeat1	gmm_cov_repeat2	allele1_num_reads	chr4_3074876_3074933_CAG_repeat_size1	chr4_3074946_3074966_CCG_repeat_size1	allele2_num_reads	chr4_3074876_3074933_CAG_repeat_size2	chr4_3074946_3074966_CCG_repeat_size2
HTT_amplicon.fastq.gz	JointGMM	2	7.8383	1.9475	705	17	10	821	55	7
```
If you copy the output to an Excel sheet, you will see the following table: 

#input_fastq | method | num_alleles | gmm_cov_repeat1 | gmm_cov_repeat2 | allele1_num_reads | chr4_3074876_3074933_CAG_repeat_size1 | chr4_3074946_3074966_CCG_repeat_size1 | allele2_num_reads | chr4_3074876_3074933_CAG_repeat_size2 | chr4_3074946_3074966_CCG_repeat_size2
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
HTT_amplicon.fastq.gz | JointGMM | 2 | 7.8383 | 1.9475 | 705 | 17 | 10 | 821 | 55 | 7


The CAG repeat sizes are 17/55. 
The CAG repeat sizes are 10/7. 
allele 1 has 705 reads
allele 2 has 821 reads 

The `*.repeat_size.txt` file reports the repeat size of each read. 

The figures are the repeat size distribution. 
	
