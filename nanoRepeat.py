#!/usr/bin/env python3

'''
Copyright (c) 2020- Children's Hospital of Philadelphia
Author: Li Fang (fangli2718@gmail.com)
              
Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

'''

import os
import sys
import time
import argparse
from argparse import RawTextHelpFormatter

import nanoRepeat_bam
import tk
from repeat_region import *

def map_fastq_to_ref_genome(in_fastq_file, ref_fasta_file, samtools, minimap2, num_cpu, bam_prefix):
    
    sam_file        = f'{bam_prefix}.sam'
    bam_file        = f'{bam_prefix}.bam'
    sorted_bam_file = f'{bam_prefix}.sorted.bam'

    if os.path.exists(sam_file):
        os.remove(sam_file)

    if os.path.exists(bam_file):
        os.remove(bam_file)

    if os.path.exists(sorted_bam_file):
        os.remove(sorted_bam_file)

    sleep_time = 5

    cmd = f'{minimap2} -ax map-ont -t {num_cpu} {ref_fasta_file} {in_fastq_file} > {sam_file} 2> /dev/null'
    tk.run_system_cmd(cmd)
    time.sleep(sleep_time)

    cmd = f'{samtools} view -hb -@ {num_cpu} {sam_file} > {bam_file} 2> /dev/null'
    tk.run_system_cmd(cmd)
    time.sleep(sleep_time)

    if os.path.getsize(bam_file) > 0:
        os.remove(sam_file)
    else:
        tk.eprint(f'ERROR! {bam_file} is empty')
        sys.exit(1)

    cmd = f'{samtools} sort -@ {num_cpu} -o {sorted_bam_file} {bam_file} 2> /dev/null'
    tk.run_system_cmd(cmd)
    time.sleep(sleep_time)

    if os.path.getsize(sorted_bam_file) > 0:
        os.remove(bam_file)
    else:
        tk.eprint(f'ERROR! {sorted_bam_file} is empty')
        sys.exit(1)

    cmd = f'{samtools} index {sorted_bam_file}'
    tk.run_system_cmd(cmd)
    time.sleep(sleep_time)

    return f'{sorted_bam_file}'

def preprocess_fastq(input_args):
    
    bam_prefix    =  f'{input_args.out_prefix}.minimap2'
    in_bam_file   = map_fastq_to_ref_genome(input_args.input, input_args.ref_fasta, input_args.samtools, input_args.minimap2, input_args.num_cpu, bam_prefix)
    
    nanoRepeat_bam.nanoRepeat_bam(input_args, in_bam_file)
    return

def main():

    program = 'nanoRepeat.py'
    examples  = f'Examples: \n'
    examples += f'\t1) python {program} -i input.bam   -t bam   -r hg38.fasta -b hg38.repeats.bed -c 4 -p prefix/of/output/files\n'
    examples += f'\t2) python {program} -i input.fastq -t fastq -r hg38.fasta -b hg38.repeats.bed -c 4 -p prefix/of/output/files\n'
    examples += f'\t3) python {program} -i input.fasta -t fasta -r hg38.fasta -b hg38.repeats.bed -c 4 -p prefix/of/output/files\n'
    
    parser = argparse.ArgumentParser(prog = program, description=f'NanoRepeat: short tandem repeat (STR) quantification from Nanopore long-read sequencing', epilog=examples, formatter_class=RawTextHelpFormatter)

    # required
    parser.add_argument('-i', '--input', required = True, metavar = 'input_file', type = str, help = '(required) path to input file (supported format: sorted_bam, fastq or fasta)')
    parser.add_argument('-t', '--type', required = True, metavar = 'input_type', type = str, help = '(required) input file type (valid values: bam, fastq or fasta)')
    parser.add_argument('-r', '--ref_fasta', required = True, metavar = 'ref.fasta',   type = str, help = '(required) path to reference genome sequence in FASTA format')
    parser.add_argument('-b', '--repeat_region_bed', required = True, metavar = 'repeat_regions.bed', type = str, help = '(required) path to the repeat region file (tab delimited text file with 4 columns: chrom start end repeat_unit_seq. Positions start from 0. Start position is self-inclusive but end position is NOT self-inclusive)')
    parser.add_argument('-o', '--out_prefix', required = True, metavar = 'path/to/out_dir/prefix_of_file_names',   type = str, help = '(required) prefix of output files')
    
    # optional 
    parser.add_argument('-c', '--num_cpu', required = False, metavar = 'INT',   type = int, default = 1,  help ='(optional) number of CPU cores (default: 1)')
    parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools',  type = str, default = 'samtools', help ='(optional) path to samtools (default: using environment default)')
    parser.add_argument('--minimap2', required = False, metavar = 'path/to/minimap2',  type = str, default = 'minimap2', help ='(optional) path to minimap2 (default: using environment default)')
    parser.add_argument('--ploidy',   required = False, metavar = 'INT', type = int, default = 2,  help ='(optional) ploidy of the sample (default: 2)')
    parser.add_argument('--anchor_len', required = False, metavar = 'INT', type = int, default = 1000, help ='(optional) length of up/downstream sequence to help identify the repeat region (default: 256 bp, increase this value if the 1000 bp up/downstream sequences are also repeat)')
    parser.add_argument('--error_rate',   required = False, metavar = 'FLOAT',  type = float, default = 0.1,  help = 'sequencing error rate (default: 0.1)')
    parser.add_argument('--max_mutual_overlap', required = False, metavar = 'FLOAT',  type = float, default = 0.15,  help = 'max mutual overlap of two alleles in terms of repeat size distribution (default value: 0.1). If the Gaussian distribution of two alleles have more overlap than this value, the two alleles will be merged into one allele.')
    parser.add_argument('--remove_noisy_reads', required = False, action='store_true', help = 'remove noisy components when there are more components than ploidy')
    parser.add_argument('--max_num_components', required = False, metavar = 'INT',  type = int, default = -1,  help = 'max number of components for the Gaussian mixture model (default value: ploidy + 20). Some noisy reads and outlier reads may form a component. Therefore the number of components is usually larger than ploidy. If your sample have too many outlier reads, you can increase this number.')

    
    if len(sys.argv) < 2 or sys.argv[1] in ['help', 'h', '-help', 'usage']:
        input_args = parser.parse_args(['--help'])
    else:
        input_args = parser.parse_args()

    input_args.type = input_args.type.lower()
    if input_args.type not in ['bam', 'fastq', 'fasta']:
        tk.eprint(f'ERROR! unknown input type: {input_args.type} valid values are: bam, fastq, fasta')
        sys.exit(1)
    
    if input_args.ploidy < 1:
        tk.eprint(f'ERROR! ploidy should be at least 1')
        sys.exit(1)

    if input_args.error_rate >= 1.0:
        tk.eprint('ERROR! --error_rate must be < 1\n')
        sys.exit(1)
    
    if input_args.max_mutual_overlap >= 1.0:
        tk.eprint('ERROR! --max_mutual_overlap must be < 1\n')
        sys.exit(1)

    if input_args.max_num_components == -1:
        input_args.max_num_components = input_args.ploidy + 20

    out_dir, out_file_prefix = os.path.split(input_args.out_prefix)
    out_dir = os.path.abspath(out_dir)
    if out_file_prefix == '':
        out_file_prefix = os.path.split(input_args.input)[1]
        input_args.out_prefix = os.path.join(out_dir, out_file_prefix)
    else:
        input_args.out_prefix = os.path.abspath(input_args.out_prefix)

    input_args.input             = os.path.abspath(input_args.input)
    input_args.ref_fasta         = os.path.abspath(input_args.ref_fasta)
    input_args.out_prefix        = os.path.abspath(input_args.out_prefix)
    input_args.repeat_region_bed = os.path.abspath(input_args.repeat_region_bed)

    tk.eprint(f'NOTICE: Input file is: {input_args.input}')
    tk.eprint(f'NOTICE: Input type is: {input_args.type}')
    tk.eprint(f'NOTICE: Referece fasta file is: {input_args.ref_fasta}')
    tk.eprint(f'NOTICE: Output prefix is: {input_args.out_prefix}')
    tk.eprint(f'NOTICE: Repeat region bed file is: {input_args.repeat_region_bed}')

    os.makedirs(out_dir, exist_ok=True)

    if input_args.type == 'fastq':
        preprocess_fastq(input_args)
    elif input_args.type == 'fasta':
        preprocess_fastq(input_args)
    elif input_args.type == 'bam':
        nanoRepeat_bam.nanoRepeat_bam(input_args, input_args.input)
    else:
        tk.eprint(f'ERROR! Unknown input type: {input_args.type} valid values are: bam, fastq, fasta')
        sys.exit(1)
    
    return


if __name__ == '__main__':
    main()
