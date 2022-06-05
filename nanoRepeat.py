#!/usr/bin/env python3

'''
Copyright (c) 2020 Children's Hospital of Philadelphia
Author: Li Fang (fangli80@foxmail.com)
              
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
import argparse
from argparse import RawTextHelpFormatter

import nanoRepeat_bam
import tk
from repeat_region import *

def map_fastq_to_ref_genome(in_fastq_file, ref_fasta_file, samtools, minimap2, num_cpu, bam_prefix):
    
    cmd = f'{minimap2} -ax map-ont -t {num_cpu} {ref_fasta_file} {in_fastq_file} > {bam_prefix}.sam'
    tk.run_system_cmd(cmd)

    cmd = f'{samtools} sort -@ {num_cpu} -o {bam_prefix}.sorted.bam {bam_prefix}.sam'
    tk.run_system_cmd(cmd)

    cmd = f'{samtools} index {bam_prefix}.sorted.bam'
    tk.run_system_cmd(cmd)

    os.remove(f'{bam_prefix}.sam')

    return f'{bam_prefix}.sorted.bam'

def preprocess_fastq(input_args):
    
    bam_prefix    =  f'{input_args.out_prefix}.minimap2'
    in_bam_file   = map_fastq_to_ref_genome(input_args.input, input_args.ref_fasta, input_args.samtools, input_args.minimap2, input_args.num_cpu, bam_prefix)
    
    nanoRepeat_bam.nanoRepeat_bam(input_args, in_bam_file)
    return

def main():

    program = 'nanoRepeat.py'
    examples  = f'Examples: \n'
    examples += f'\t1) python {program} -i input.bam   -t bam   -r hg38.fasta -b hg38.repeats.bed -c 4 -o prefix/of/output/files\n'
    examples += f'\t2) python {program} -i input.fastq -t fastq -r hg38.fasta -b hg38.repeats.bed -c 4 -o prefix/of/output/files\n'
    examples += f'\t3) python {program} -i input.fasta -t fasta -r hg38.fasta -b hg38.repeats.bed -c 4 -o prefix/of/output/files\n'
    
    parser = argparse.ArgumentParser(prog = program, description=f'NanoRepeat: short tandem repeat (STR) quantification from Nanopore long-read sequencing', epilog=examples, formatter_class=RawTextHelpFormatter)

    # required
    parser.add_argument('-i', '--input', required = True, metavar = 'input_file', type = str, help = '(required) path to input file (supported format: sorted_bam, fastq or fasta)')
    parser.add_argument('-t', '--type', required = True, metavar = 'input_type', type = str, help = '(required) input file type (valid values: bam, fastq or fasta)')
    parser.add_argument('-r', '--ref_fasta', required = True, metavar = 'ref.fasta',   type = str, help = '(required) path to reference genome sequence in FASTA format')
    parser.add_argument('-b', '--repeat_region_bed', required = True, metavar = 'repeat_regions.bed', type = str, help = '(required) path to repeat region file (in bed format)')
    parser.add_argument('-o', '--out_prefix', required = True, metavar = 'prefix/of/output/files',   type = str, help = '(required) prefix of output files')
    
    # optional 
    parser.add_argument('-c', '--num_cpu', required = False, metavar = 'INT',   type = int, default = 1,  help ='(optional) number of CPU cores (default: 1)')
    parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools',  type = str, default = 'samtools', help ='(optional) path to samtools (default: using environment default)')
    parser.add_argument('--minimap2', required = False, metavar = 'path/to/minimap2',  type = str, default = 'minimap2', help ='(optional) path to minimap2 (default: using environment default)')
    parser.add_argument('--ploidy',   required = False, metavar = 'INT', type = int, default = 2,  help ='(optional) ploidy of the sample (default: 2)')
    parser.add_argument('--anchor_len', required = False, metavar = 'INT', type = int, default = 1000, help ='(optional) length of up/downstream sequence to help identify the repeat region (default: 256 bp, increase this value if the 1000 bp up/downstream sequences are also repeat)')

    
    if len(sys.argv) < 2 or sys.argv[1] in ['help', 'h', '-help', 'usage']:
        input_args = parser.parse_args(['--help'])
    else:
        input_args = parser.parse_args()

    input_args.type = input_args.type.lower()
    if input_args.type not in ['bam', 'fastq', 'fasta']:
        tk.eprint(f'ERROR! unknown input type: {input_args.type} valid values are: bam, fastq, fasta')
        sys.exit(1)
    
    input_args.input             = os.path.abspath(input_args.input)
    input_args.ref_fasta         = os.path.abspath(input_args.ref_fasta)
    input_args.out_prefix        = os.path.abspath(input_args.out_prefix)
    input_args.repeat_region_bed = os.path.abspath(input_args.repeat_region_bed)

    tk.eprint(f'NOTICE: input file is: {input_args.input}')
    tk.eprint(f'NOTICE: input type is: {input_args.type}')
    tk.eprint(f'NOTICE: referece fasta file is: {input_args.ref_fasta}')
    tk.eprint(f'NOTICE: output prefix is: {input_args.out_prefix}')
    tk.eprint(f'NOTICE: repeat region bed file is: {input_args.repeat_region_bed}')

    out_dir = os.path.split(input_args.out_prefix)[0]
    os.makedirs(out_dir, exist_ok=True)

    if input_args.type == 'fastq':
        preprocess_fastq(input_args)
    elif input_args.type == 'fasta':
        preprocess_fastq(input_args)
    elif input_args.type == 'bam':
        nanoRepeat_bam.nanoRepeat_bam(input_args, input_args.input)
    else:
        tk.eprint(f'ERROR! unknown input type: {input_args.type} valid values are: bam, fastq, fasta')
        sys.exit(1)
    
    return


if __name__ == '__main__':
    main()
