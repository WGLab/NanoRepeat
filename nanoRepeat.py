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
import nanoRepeat_fq

def main():

    program = 'nanoRepeat.py'
    examples  = f'Examples: \n'
    examples += f'\t1) bam input:   python {program} bam   -i input.bam   -f hg19.fasta -r repeats.txt -o ./output\n'
    examples += f'\t2) fastq input: python {program} fastq -i input.fastq -f hg19.fasta -r repeats.txt -o ./output\n'

    parser = argparse.ArgumentParser(prog = program, description=f'A toolkit for short tandem repeat (STR) quantification from Nanopore long-read sequencing data.', epilog=examples, formatter_class=RawTextHelpFormatter)
    
    subparsers = parser.add_subparsers(dest="subparser_name") # this line changed
    bam_parser = subparsers.add_parser('bam')
    bam_parser.add_argument('-i', '--in_bam', required = True, metavar = 'input.bam',   type = str, help = '(required) path to input bam file (must be sorted by coordinates)')
    bam_parser.add_argument('-f', '--ref_fasta', required = True, metavar = 'ref.fasta',   type = str, help = '(required) path to reference genome sequence in FASTA format')
    bam_parser.add_argument('-r', '--repeat_regions', required = True, metavar = 'repeat_regions.bed', type = str, help = '(required) path to repeat region file (in bed format)')
    bam_parser.add_argument('-o', '--out_dir', required = True, metavar = 'path/to/out_dir',   type = str, help = '(required) path to the output directory')
    bam_parser.add_argument('-t', '--num_threads', required = False, metavar = 'INT',   type = int, default = 1,  help ='(optional) number of threads used by minimap2 (default: 1)')
    bam_parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools',  type = str, default = 'samtools', help ='(optional) path to samtools (default: using environment default)')
    bam_parser.add_argument('--minimap2', required = False, metavar = 'path/to/minimap2',  type = str, default = 'minimap2', help ='(optional) path to minimap2 (default: using environment default)')
    bam_parser.add_argument('--ploidy',   required = False, metavar = 'INT',   type = int, default = 2,  help ='(optional) ploidy of the sample (default: 2)')
    bam_parser.add_argument('--anchor_len', required = False, metavar = 'INT',   type = int, default = 256, help ='(optional) length of up/downstream sequence to help identify the repeat region (default: 256 bp, increase this value if the 1000 bp up/downstream sequences are also repeat)')

    fastq_parser = subparsers.add_parser('fastq')
    fastq_parser.add_argument('-i', '--in_fastq', required = True, metavar = 'input.fastq',   type = str, help = 'path to input fastq file')
    fastq_parser.add_argument('-f', '--ref_fasta', required = True, metavar = 'ref.fasta',   type = str, help = '(required) path to reference genome sequence in FASTA format')
    fastq_parser.add_argument('-r', '--repeat_regions', required = True, metavar = 'repeat_regions.txt', type = str, help = '(required) path to repeat region file')
    fastq_parser.add_argument('-o', '--out_dir', required = True, metavar = 'path/to/out_dir',   type = str, help = '(required) path to the output directory')
    fastq_parser.add_argument('-t', '--num_threads', required = False, metavar = 'INT',   type = int, default = 1,  help ='(optional) number of threads used by minimap2 (default: 1)')
    fastq_parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools',  type = str, default = 'samtools', help ='(optional) path to samtools (default: using environment default)')
    fastq_parser.add_argument('--minimap2', required = False, metavar = 'path/to/minimap2',  type = str, default = 'minimap2', help ='(optional) path to minimap2 (default: using environment default)')
    fastq_parser.add_argument('--ploidy',   required = False, metavar = 'INT',   type = int, default = 2,  help ='(optional) ploidy of the sample (default: 2)')
    fastq_parser.add_argument('--anchor_len', required = False, metavar = 'INT',   type = int, default = 256, help ='(optional) length of up/downstream sequence to help identify the repeat region (default: 256 bp, increase this value if the 1000 bp up/downstream sequences are also repeat)')
    
    if len(sys.argv) < 2 or sys.argv[1] in ['help', 'h', '-h', '-help', '--help', 'usage']:
        input_args = parser.parse_args(['--help'])
    else:
        input_args = parser.parse_args()

    if input_args.subparser_name == 'bam':
        nanoRepeat_bam.nanoRepeat_bam(input_args)
    elif input_args.subparser_name == 'fastq':
        nanoRepeat_fq.nanoRepeat_fq(input_args)
    else:
        input_args = parser.parse_args(['--help'])

    return


    

if __name__ == '__main__':
    main()
