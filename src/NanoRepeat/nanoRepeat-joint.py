#!/usr/bin/env python3


'''
Copyright (c) 2020-2023 Children's Hospital of Philadelphia
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
from typing import final
import numpy as np

import argparse
import shutil


from NanoRepeat import tk
from NanoRepeat.paf import *
from NanoRepeat.split_alleles import *

class Repeat:
    def __init__(self):
        self.repeat_id = ''
        self.chrom = ''
        self.start = -1
        self.end = -1
        self.repeat_unit = ''
        self.repeat_unit_size = 0
        self.min_size = 0
        self.max_size = 1000
        self.round1_min_size = 0
        self.round1_max_size = 1000
        self.is_main = 0
    
    def init_from_string(self, string):
        col_list = string.split(':')
        if len(col_list) != 5:
            tk.eprint('ERROR! --repeat1 and --repeat2 should be in this format: chr:start:end:repeat_unit:max_size (coordinates are 0-based, e.g. chr4:3074876:3074933:CAG:200).')
            sys.exit(1)
        
        self.chrom, self.start, self.end, self.repeat_unit, self.max_size = col_list
        self.start = int(self.start)
        self.end   = int(self.end)
        self.repeat_unit_size = len(self.repeat_unit)
        self.min_size = 0
        self.max_size = int(self.max_size)
        self.repeat_id = '-'.join(col_list[0:4])
        return self

class Round1Estimation:
    def __init__(self):
        self.repeat1_count_range_dict = dict()
        self.repeat2_count_range_dict = dict()
        self.potential_repeat_region_dict = dict()
        self.bad_reads_set = set()
    
    def output_repeat2_boundaries(self):
        out_string = ''
        for readname in self.repeat2_count_range_dict:
            lower_bound, upper_bound = self.repeat2_count_range_dict[readname]
            out_string += f'{readname}\t{lower_bound}\t{upper_bound}\n'

        return out_string

class RepeatSize:
    def __init__(self):
        self.repeat1_count_dict = dict()
        self.repeat2_count_dict = dict()
        self.step_size1 = 1
        self.step_size2 = 1

def parse_user_arguments():

    parser = argparse.ArgumentParser(description='Joint quantification of two adjacent tandem repeats from long-read amplicon sequencing data ')
    ### required arguments ###
    parser.add_argument('-i', '--in_fq',        required = True,  metavar = 'PATH',   type = str, help = 'input fastq file')
    parser.add_argument('-r', '--ref_fasta',    required = True,  metavar = 'PATH',   type = str, help = 'reference genome sequence in FASTA format')
    parser.add_argument('-1', '--repeat1',      required = True,  metavar = 'chr:start:end:repeat_unit:max_size', type = str, help = 'first tandem repeat in the format of chr:start:end:repeat_unit:max_size. Positions start from 0. Start position is self-inclusive but end position is NOT self-inclusive (e.g. chr4:3074876:3074933:CAG:200)')
    parser.add_argument('-2', '--repeat2',      required = True,  metavar = 'chr:start:end:repeat_unit:max_size', type = str, help = 'second tandem repeat in the format of chr:start:end:repeat_unit:max_size. Positions start from 0. Start position is self-inclusive but end position is NOT self-inclusive (e.g. chr4:3074946:3074966:CCG:20)')
    parser.add_argument('-o', '--out_prefix',   required = True, metavar = 'path/to/out_dir/prefix_of_file_names',   type = str, help = '(required) prefix of output files.')
    

    ### optional arguments ### ploidy
    parser.add_argument('-d', '--data_type', required = False, metavar = 'data_type',  type = str, default = 'ont', help = '(required) sequencing data type. Should be one of the following: ont, ont_sup, ont_q20, clr, hifi')
    parser.add_argument('-c', '--num_threads',  required = False, metavar = 'INT',    type = int, default = 1,  help = 'number of threads used by minimap2 (default: 1)')
    parser.add_argument('--minimap2',     required = False, metavar = 'PATH',   type = str, default = '', help = 'path to minimap2 (default: using environment default)')
    parser.add_argument('--ploidy',       required = False, metavar = 'INT',    type = int, default = 2,  help = 'ploidy of the sample (default: 2)')
    parser.add_argument('--error_rate',   required = False, metavar = 'FLOAT',  type = float, default = 0.1,  help = 'sequencing error rate (default: 0.1)')
    parser.add_argument('--max_mutual_overlap', required = False, metavar = 'FLOAT',  type = float, default = 0.1,  help = 'max mutual overlap of two alleles in terms of repeat size distribution (default value: 0.1). If the Gaussian distribution of two alleles have more overlap than this value, the two alleles will be merged into one allele.')
    parser.add_argument('--remove_noisy_reads', required = False, action='store_true', help = 'remove noisy components when there are more components than ploidy')
    parser.add_argument('--max_num_components', required = False, metavar = 'INT',  type = int, default = -1,  help = 'max number of components for the Gaussian mixture model (default value: ploidy + 20). Some noisy reads and outlier reads may form a component. Therefore the number of components is usually larger than ploidy. If your sample have too many outlier reads, you can increase this number.')


    if len(sys.argv) < 2 or sys.argv[1] in ['help', 'h', '-help', 'usage']:
        input_args = parser.parse_args(['--help'])
    else:
        input_args = parser.parse_args()

    if input_args.data_type not in ['ont', 'ont_sup', 'ont_q20', 'clr', 'hifi']:
        tk.eprint(f'ERROR! data_type should be one of the following: ont, ont_sup, clr, hifi\n')
        tk.eprint(f'ont:     Oxford Nanopore sequencing, NOT Super Accuracy mode\n')
        tk.eprint(f'ont_sup: Oxford Nanopore sequencing, Super Accuracy mode\n')
        tk.eprint(f'ont_q20: Oxford Nanopore sequencing, Q20+ chemistry\n')
        tk.eprint(f'clr:     PacBio sequencing, Continuous Long Reads (CLR)\n')
        tk.eprint(f'hifi:    PacBio sequencing, HiFi/CCS reads\n')
        sys.exit(1)
        
    if input_args.ploidy < 1:
        tk.eprint('ERROR: --ploidy must be >= 1 !\n')
        sys.exit(1)

    if input_args.error_rate >= 1.0:
        tk.eprint('ERROR! --error_rate must be < 1\n')
        sys.exit(1)
    
    if input_args.max_mutual_overlap >= 1.0:
        tk.eprint('ERROR! --max_mutual_overlap must be < 1\n')
        sys.exit(1)

    if input_args.max_num_components == -1:
        input_args.max_num_components = input_args.ploidy + 20
    
    tk.check_input_file_exists(input_args.in_fq)
    tk.check_input_file_exists(input_args.ref_fasta)

    input_args.in_fq       = os.path.abspath(input_args.in_fq)
    input_args.ref_fasta   = os.path.abspath(input_args.ref_fasta)

    if input_args.minimap2 == '': 
        input_args.minimap2 = tk.find_executable_path('minimap2')
        if not input_args.minimap2:
            tk.eprint('ERROR! minimap2 was not found! Please supply the path to --minimap2')
            sys.exit(1)
        else:
            if os.path.exists(input_args.minimap2):
                tk.eprint('NOTICE: Found path to minimap2: %s' % input_args.minimap2)
            else:
                tk.eprint('ERROR! minimap2 was not found! Please supply the path to --minimap2')

    tk.check_input_file_exists(input_args.minimap2)

    if len(input_args.repeat1.split(':')) != 5:
        tk.eprint('ERROR! --repeat1 should be in this format: chr:start:end:repeat_unit:max_size (coordinates are 0-based, e.g. chr4:3074876:3074933:CAG:200)')
        sys.exit(1)
    
    if len(input_args.repeat2.split(':')) != 5:
        tk.eprint('ERROR! --repeat2 should be in this format: chr:start:end:repeat_unit:max_size (coordinates are 0-based, e.g. chr4:3074946:3074966:CCG:20)')
        sys.exit(1)

    return input_args

def main():

    input_args = parse_user_arguments()
    nanoRepeat_joint (input_args)

    return

def nanoRepeat_joint (input_args):

    n_input_reads = tk.count_fastq(input_args.in_fq)
    if n_input_reads < input_args.ploidy:
        tk.eprint(f'ERROR: No enough reads for analysis. Ploidy was set to {input_args.ploidy} but there were only {n_input_reads} reads in the fastq file: {input_args.in_fq}\n')
        sys.exit(1)
    
    out_dir, out_file_prefix = os.path.split(input_args.out_prefix)
    out_dir = os.path.abspath(out_dir)

    if out_file_prefix == '':
        out_file_prefix = os.path.split(input_args.in_fq)[1]
        input_args.out_prefix = os.path.join(out_dir, out_file_prefix)
    else:
        input_args.out_prefix = os.path.abspath(input_args.out_prefix)

    tk.create_dir(out_dir)

    tk.eprint(f'NOTICE: Input file is: {input_args.in_fq }')
    tk.eprint(f'NOTICE: Referece fasta file is: {input_args.ref_fasta}')
    tk.eprint(f'NOTICE: Repeat1 is: {input_args.repeat1}')
    tk.eprint(f'NOTICE: Repeat2 is: {input_args.repeat2}')
    tk.eprint(f'NOTICE: Output prefix is: {input_args.out_prefix}')
    

    repeat1 = Repeat().init_from_string(input_args.repeat1)
    repeat2 = Repeat().init_from_string(input_args.repeat2)

    if repeat1.chrom != repeat2.chrom:
        tk.eprint('ERROR: joint quantification only works with two nearby repeat regions. The two repeats in the config_file are in different chromosomes and there is no need to do joint quantification.\n')
        sys.exit(1)

    if repeat1.start > repeat2.start:
        repeat1, repeat2 = tk.switch_two_objects(repeat1, repeat2)

    repeat1.max_size += 10
    repeat2.max_size += 10
    max_anchor_len = 1000
    
    if repeat1.end + 100 < repeat2.start:
        tk.eprint('ERROR: joint quantification only works with two nearby repeat regions (distance < 100 bp). The two repeats are far away from each other and there is no need to do joint quantification.\n')
        sys.exit(1)

    temp_out_dir = os.path.join(out_dir, f'{out_file_prefix}.NanoRepeat.temp')
    tk.create_dir(temp_out_dir)

    repeat_chrom_seq = tk.read_one_chr_from_fasta_file(input_args.ref_fasta, repeat1.chrom)
    if len(repeat_chrom_seq) == 0:
        tk.eprint('ERROR: ref_fasta file: %s has no valid sequence!\n' % input_args.ref_fasta)
        sys.exit(1)

    initial_estimation = initial_estimate_repeat_size(input_args.minimap2, repeat_chrom_seq, input_args.in_fq, input_args.data_type, input_args.num_threads, repeat1, repeat2, max_anchor_len, temp_out_dir)

    final_estimation = fine_tune_read_count(initial_estimation, input_args.in_fq, repeat_chrom_seq, repeat1, repeat2, input_args.minimap2, input_args.data_type, input_args.num_threads, temp_out_dir)
    
    read_repeat_joint_count_dict = output_repeat_size_2d (input_args.in_fq, repeat1.repeat_id, repeat2.repeat_id, input_args.out_prefix, final_estimation.repeat1_count_dict, final_estimation.repeat2_count_dict)

    tk.eprint('NOTICE: Phasing reads using GMM')
    num_removed_reads = 0
    split_alleles_using_gmm_2d (input_args.ploidy, input_args.error_rate, input_args.max_mutual_overlap, input_args.remove_noisy_reads, input_args.max_num_components, repeat1, repeat2, read_repeat_joint_count_dict, num_removed_reads, input_args.in_fq, input_args.out_prefix)

    shutil.rmtree(temp_out_dir)
    tk.eprint(f'NOTICE: Program finished. Output files are here: {out_dir}\n')

    return


def fine_tune_read_count(initial_estimation, in_fastq_file, repeat_chrom_seq, repeat1, repeat2, minimap2, data_type, num_threads, out_dir):

    assert repeat1.chrom == repeat2.chrom
    assert repeat1.start < repeat2.start

    repeat1.round1_max_size = 0
    repeat2.round1_max_size = 0
    repeat1.round1_min_size = repeat1.max_size
    repeat2.round1_min_size = repeat2.max_size
    for readname, value in initial_estimation.repeat1_count_range_dict.items():
        min_repeat_size, max_repeat_size = value
        if max_repeat_size > repeat1.round1_max_size:
            repeat1.round1_max_size = max_repeat_size
        if min_repeat_size < repeat1.round1_min_size:
            repeat1.round1_min_size = min_repeat_size

    for readname, value in initial_estimation.repeat2_count_range_dict.items():
        min_repeat_size, max_repeat_size = value
        if max_repeat_size > repeat2.round1_max_size:
            repeat2.round1_max_size = max_repeat_size
        if min_repeat_size < repeat2.round1_min_size:
            repeat2.round1_min_size = min_repeat_size

    
    if repeat1.round1_max_size > repeat1.max_size: repeat1.round1_max_size = repeat1.max_size
    if repeat2.round1_max_size > repeat2.max_size: repeat2.round1_max_size = repeat2.max_size

    tk.eprint('NOTICE: In round 1 estimation, repeat 1 (%s) is in the range of (%d, %d)' % (repeat1.repeat_unit, repeat1.round1_min_size, repeat1.round1_max_size))
    tk.eprint('NOTICE: In round 1 estimation, repeat 2 (%s) is in the range of (%d, %d)' % (repeat2.repeat_unit, repeat2.round1_min_size, repeat2.round1_max_size))
    
    fastq_dict = fastq_file_to_dict(in_fastq_file)

    round2_estimation = round2_estimation_of_repeat_size(initial_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, data_type, num_threads, out_dir)
    
    if round2_estimation.step_size1 > 1 and round2_estimation.step_size2 > 1:
        final_estimation = round3_estimation_of_repeat_size(initial_estimation, round2_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, data_type, num_threads, out_dir)
    else:
        final_estimation = round2_estimation
    
    return final_estimation

def round3_estimation_of_repeat_size(initial_estimation, round2_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, data_type, num_threads, out_dir):

    if len(round2_estimation.repeat1_count_dict) == 0 or len(round2_estimation.repeat2_count_dict) == 0:
        return RepeatSize()
    
    assert repeat1.chrom == repeat2.chrom
    assert repeat1.start < repeat2.start
    max_flanking_len = 1000

    buffer_size1 = round2_estimation.step_size1
    buffer_size2 = round2_estimation.step_size2
    tk.eprint(f'NOTICE: Fine-tuning repeat size. buffer_size1 = {buffer_size1}, buffer_size2 = {buffer_size2}')

    round2_size1_list = list()
    round2_size2_list = list()

    for readname in round2_estimation.repeat1_count_dict:
        if readname not in round2_estimation.repeat2_count_dict: continue
        size1 = round2_estimation.repeat1_count_dict[readname]
        size2 = round2_estimation.repeat2_count_dict[readname]
        round2_size1_list.append(size1)
        round2_size2_list.append(size2)

    min_size1 = int(min(round2_size1_list) - buffer_size1)
    max_size1 = int(max(round2_size1_list) + buffer_size1 + 2)
    if min_size1 < 0: min_size1 = 0
    min_size2 = int(min(round2_size2_list) - buffer_size2)
    max_size2 = int(max(round2_size2_list) + buffer_size2 + 2)
    if min_size2 < 0: min_size2 = 0

    ## round 3 alignment ##
    left_anchor_seq, mid_anchor_seq, right_anchor_seq = extract_anchor_seq_for_two_repeats (repeat_chrom_seq, repeat1, repeat2, max_flanking_len)
    left_anchor_len  = len(left_anchor_seq)
    mid_anchor_len   = len(mid_anchor_seq)

    preset = tk.get_preset_for_minimap2(data_type)
    round3_paf_file =  os.path.join(out_dir, 'round3.paf')
    round3_paf_f = open(round3_paf_file, 'w')
    round3_paf_f.close()

    for repeat_count1 in range(min_size1, max_size1):
        for repeat_count2 in range(min_size2, max_size2):
            tmp_fastq_file = os.path.join(out_dir, '%d.%d.round3.fastq' % (repeat_count1, repeat_count2) )
            tmp_fastq_f = open(tmp_fastq_file, 'w')
            n_tmp_reads = 0
            for readname in fastq_dict:
                if readname not in round2_estimation.repeat1_count_dict: continue
                if readname not in round2_estimation.repeat2_count_dict: continue
                size1 = round2_estimation.repeat1_count_dict[readname]
                size2 = round2_estimation.repeat2_count_dict[readname]
                if repeat_count1 < size1 - buffer_size1 or repeat_count1 >= size1 + buffer_size1: continue
                if repeat_count2 < size2 - buffer_size2 or repeat_count2 >= size2 + buffer_size2: continue
                r1min1, r1max1 = initial_estimation.repeat1_count_range_dict[readname]
                r1min2, r1max2 = initial_estimation.repeat2_count_range_dict[readname]
                if repeat_count1 < r1min1 or repeat_count1 >= r1max1: continue
                if repeat_count2 < r1min2 or repeat_count2 >= r1max2: continue

                tmp_fastq_f.write(fastq_dict[readname])
                n_tmp_reads += 1
            tmp_fastq_f.close()
            if n_tmp_reads == 0: 
                tk.rm(tmp_fastq_file)
                continue
            tmp_ref_file = os.path.join(out_dir, '%d.%d.ref.fasta' % (repeat_count1, repeat_count2))
            build_fasta_template_for_two_repeats(left_anchor_seq, mid_anchor_seq, right_anchor_seq, repeat1, repeat2, repeat_count1, repeat_count2, tmp_ref_file)
            cmd = f'{minimap2} -c --eqx -t {num_threads} {preset} {tmp_ref_file} {tmp_fastq_file} >> {round3_paf_file} 2> /dev/null'
            tk.run_system_cmd(cmd)
            tk.rm(tmp_ref_file)
            tk.rm(tmp_fastq_file)

    round3_estimation = estimate_two_repeats_from_paf(round3_paf_file, left_anchor_len, mid_anchor_len, repeat1, repeat2)
    round3_estimation.step_size1 = 1
    round3_estimation.step_size2 = 1

    return round3_estimation

def choose_best_step_size(repeat, count_range_dict):

    max_len = 50
    max_step_size = int(max_len/repeat.repeat_unit_size)
    if max_step_size < 1: max_step_size = 1

    error_list = list()
    for readname in count_range_dict:
        a, b = count_range_dict[readname]
        error = b - a
        error_list.append(error)

    l = np.mean(error_list)

    count_list = list()
    for size in range(1, max_step_size+1):
        count = int(l/size)+1
        count += size * 2 + 2
        tup = (size, count)
        count_list.append(tup)

    count_list.sort(key = lambda x:x[1])
  
    return count_list[0][0]

def round2_estimation_of_repeat_size(initial_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, data_type, num_threads, out_dir):
    
    assert repeat1.chrom == repeat2.chrom
    assert repeat1.start < repeat2.start
    max_flanking_len = 1000

    step_size1 = choose_best_step_size(repeat1, initial_estimation.repeat1_count_range_dict)
    step_size2 = choose_best_step_size(repeat2, initial_estimation.repeat2_count_range_dict)
    
    tk.eprint(f'NOTICE: Round 2 estimation. step_size1 = {step_size1}; step_size2 = {step_size2}')

    ## round 2 alignment ##
    left_anchor_seq, mid_anchor_seq, right_anchor_seq = extract_anchor_seq_for_two_repeats (repeat_chrom_seq, repeat1, repeat2, max_flanking_len)
    left_anchor_len  = len(left_anchor_seq)
    mid_anchor_len   = len(mid_anchor_seq)

    preset = tk.get_preset_for_minimap2(data_type)
    round2_paf_file =  os.path.join(out_dir, 'round2.paf')
    round2_paf_f = open(round2_paf_file, 'w')
    round2_paf_f.close()

    for repeat_count1 in range(repeat1.round1_min_size, repeat1.round1_max_size + 1, step_size1):
        for repeat_count2 in range(repeat2.round1_min_size, repeat2.round1_max_size + 1, step_size2):
            tmp_fastq_file = os.path.join(out_dir, '%d.%d.round2.fastq' % (repeat_count1, repeat_count2) )
            tmp_fastq_f = open(tmp_fastq_file, 'w')
            n_tmp_reads = 0
            for readname in fastq_dict:
                if readname not in initial_estimation.repeat1_count_range_dict: continue
                if readname not in initial_estimation.repeat2_count_range_dict: continue
                min1, max1 = initial_estimation.repeat1_count_range_dict[readname]
                min2, max2 = initial_estimation.repeat2_count_range_dict[readname]
                if repeat_count1 >= min1 and repeat_count1 < max1 and repeat_count2 >= min2 and repeat_count2 < max2:                
                    tmp_fastq_f.write(fastq_dict[readname])
                    n_tmp_reads += 1
            tmp_fastq_f.close()
            if n_tmp_reads == 0:
                tk.rm(tmp_fastq_file)
                continue
            tmp_ref_file = os.path.join(out_dir, '%d.%d.ref.fasta' % (repeat_count1, repeat_count2))
            build_fasta_template_for_two_repeats(left_anchor_seq, mid_anchor_seq, right_anchor_seq, repeat1, repeat2, repeat_count1, repeat_count2, tmp_ref_file)
            cmd = f'{minimap2} {preset} -c --eqx -t {num_threads}  {tmp_ref_file} {tmp_fastq_file} >> {round2_paf_file} 2> /dev/null'
            tk.run_system_cmd(cmd)
            tk.rm(tmp_ref_file)
            tk.rm(tmp_fastq_file)

    round2_estimation = estimate_two_repeats_from_paf(round2_paf_file, left_anchor_len, mid_anchor_len, repeat1, repeat2)
    round2_estimation.step_size1 = step_size1
    round2_estimation.step_size2 = step_size2

    return round2_estimation

def estimate_two_repeats_from_paf(in_paf_file, left_anchor_len, mid_anchor_len, repeat1, repeat2):

    repeat_size_estimation = RepeatSize()

    in_paf_f = open(in_paf_file, 'r')
    read_align_score_dict = dict()
    while 1:
        line = in_paf_f.readline()
        if not line: break
        line = line.strip()
        if not line: continue

        col_list = line.strip().split('\t')
        paf = PAF(col_list)
        repeat_size1, repeat_size2 = paf.tname.split('-')
        repeat_size1 = int(repeat_size1)
        repeat_size2 = int(repeat_size2)

        a = left_anchor_len - 10
        if a < 0: a = 0
        b = left_anchor_len + repeat1.repeat_unit_size * repeat_size1 + mid_anchor_len + repeat2.repeat_unit_size * repeat_size2 + 10
        if b > paf.tlen: b = paf.tlen
        repeat_region_align_score = tk.target_region_alignment_stats_from_cigar(paf.cigar, paf.tstart, paf.tend, a, b).score
    
        tup = (repeat_size1, repeat_size2, repeat_region_align_score, paf.align_score, paf.qlen)
        if paf.qname not in read_align_score_dict:
            read_align_score_dict[paf.qname] = list()
        read_align_score_dict[paf.qname].append(tup)    

    in_paf_f.close()

    for readname in read_align_score_dict:
        tup_list = read_align_score_dict[readname]
        tup_list.sort(key = lambda x:x[2], reverse = True)
        max_score = tup_list[0][2]
        max_score_size1_list = list()
        max_score_size2_list = list()
        for tup in tup_list:
            if tup[2] == max_score:
                size1 = tup[0]
                size2 = tup[1]
                max_score_size1_list.append(size1)
                max_score_size2_list.append(size2)
            else:
                break
        
        round2_repeat_size1 = np.mean(max_score_size1_list)
        round2_repeat_size2 = np.mean(max_score_size2_list)
        repeat_size_estimation.repeat1_count_dict[readname] = round2_repeat_size1
        repeat_size_estimation.repeat2_count_dict[readname] = round2_repeat_size2
    
    return repeat_size_estimation

def extract_anchor_seq_for_two_repeats (repeat_chrom_seq, repeat1, repeat2, max_flanking_len):

    left_end_pos = repeat1.start
    left_start_pos = left_end_pos - max_flanking_len
    if left_start_pos < 0: left_start_pos = 0
    
    mid_start_pos = repeat1.end
    mid_end_pos = repeat2.start

    right_start_pos = repeat2.end
    right_end_pos = right_start_pos + max_flanking_len
    if right_end_pos > len(repeat_chrom_seq): right_end_pos = len(repeat_chrom_seq)

    left_anchor_seq  = repeat_chrom_seq[left_start_pos:left_end_pos]
    mid_anchor_seq   = repeat_chrom_seq[mid_start_pos:mid_end_pos]
    right_anchor_seq = repeat_chrom_seq[right_start_pos:right_end_pos]

    return left_anchor_seq, mid_anchor_seq, right_anchor_seq

def build_fasta_template_for_two_repeats(left_anchor_seq, mid_anchor_seq, right_anchor_seq, repeat1, repeat2, repeat_count1, repeat_count2, tmp_ref_file):

    tmp_ref_f = open(tmp_ref_file, 'w')
    tmp_ref_f.write('>%d-%d\n' % (repeat_count1, repeat_count2))
    seq = left_anchor_seq + repeat1.repeat_unit * repeat_count1 + mid_anchor_seq + repeat2.repeat_unit * repeat_count2 + right_anchor_seq + '\n'
    tmp_ref_f.write(seq)
    tmp_ref_f.close()

    return

def initial_estimate_repeat_size(minimap2, repeat_chrom_seq, in_fastq_file, data_type, num_threads, repeat1, repeat2, max_anchor_len, out_dir):

    tk.eprint('NOTICE: Round 1 estimation')
    assert repeat1.chrom == repeat2.chrom
    assert repeat1.start < repeat2.start

    left_anchor_seq, right_anchor_seq = tk.extract_anchor_sequence(repeat_chrom_seq, repeat1.start, repeat2.end, max_anchor_len)
    left_anchor_len  = len(left_anchor_seq)
    right_anchor_len = len(right_anchor_seq)

    ## make template fasta file ##
    left_template_fasta_file = os.path.join(out_dir, 'round1_templates.left_anchor.fasta')
    left_template_name  = 'left_anchor_%d_%d_%s' % (len(left_anchor_seq), repeat1.max_size, repeat1.repeat_unit)
    left_template_seq   = left_anchor_seq + repeat1.repeat_unit * repeat1.max_size

    left_template_fasta_f = open(left_template_fasta_file, 'w')
    left_template_fasta_f.write('>%s\n' % left_template_name)
    left_template_fasta_f.write('%s\n' % left_template_seq)
    left_template_fasta_f.close()


    right_template_fasta_file = os.path.join(out_dir, 'round1_templates.right_anchor.fasta')
    right_template_name = 'right_anchor_%d_%d_%s_revc' % (len(right_anchor_seq), repeat2.max_size, repeat2.repeat_unit)
    right_template_seq  = repeat2.repeat_unit * repeat2.max_size + right_anchor_seq
    right_template_seq  = tk.rev_comp(right_template_seq)
    
    right_template_fasta_f = open(right_template_fasta_file, 'w')
    right_template_fasta_f.write('>%s\n' % right_template_name)
    right_template_fasta_f.write('%s\n' % right_template_seq)
    right_template_fasta_f.close()

    round1_paf_file = os.path.join(out_dir, 'round1.paf')
    round1_left_anchor_paf_file  = os.path.join(out_dir, 'round1.left.paf')
    round1_right_anchor_paf_file = os.path.join(out_dir, 'round1.right.paf')
    preset = tk.get_preset_for_minimap2(data_type)

    cmd = f'{minimap2} {preset} -c --eqx -t {num_threads} {left_template_fasta_file} {in_fastq_file} > {round1_left_anchor_paf_file} 2> /dev/null'
    tk.run_system_cmd(cmd)
    
    cmd = f'{minimap2} {preset} -c --eqx -t {num_threads} {right_template_fasta_file} {in_fastq_file} > {round1_right_anchor_paf_file} 2> /dev/null'
    tk.run_system_cmd(cmd)

    cmd = f'cat {round1_left_anchor_paf_file} {round1_right_anchor_paf_file} | sort -k1 - > {round1_paf_file}'
    tk.run_system_cmd(cmd)

    tk.rm(round1_left_anchor_paf_file)
    tk.rm(round1_right_anchor_paf_file)
    initial_estimation = round1_estimation_from_paf(round1_paf_file, repeat1, repeat2, left_anchor_len, right_anchor_len)

    return initial_estimation

def round1_estimation_from_paf(round1_paf_file, repeat1, repeat2, left_anchor_len, right_anchor_len):

    initial_estimation = Round1Estimation()
    round1_paf_f = open(round1_paf_file, 'r')
    read_paf_list = list()
    while 1:
        line = round1_paf_f.readline()
        if not line: break
        line = line.strip()
        if not line: continue

        col_list = line.strip().split('\t')
        paf = PAF(col_list)

        if len(read_paf_list) == 0 or paf.qname == read_paf_list[0].qname:
            read_paf_list.append(paf)
        else:
            round1_estimation_for1read(read_paf_list, repeat1, repeat2, left_anchor_len, right_anchor_len, initial_estimation)
            read_paf_list.clear()
            read_paf_list.append(paf)

    round1_paf_f.close()
    
    round1_estimation_for1read(read_paf_list, repeat1, repeat2, left_anchor_len, right_anchor_len, initial_estimation)

    for readname in initial_estimation.bad_reads_set:
        initial_estimation.repeat1_count_range_dict.pop(readname, None)
        initial_estimation.repeat2_count_range_dict.pop(readname, None)

    return initial_estimation
            

def round1_estimation_for1read(read_paf_list, repeat1, repeat2, left_anchor_len, right_anchor_len, initial_estimation):

    left_boundary_pos  = left_anchor_len
    right_boundary_pos = right_anchor_len
    
    left_anchor_paf_list  = list()
    right_anchor_paf_list = list()
    
    for paf in read_paf_list:
        if paf.mapq < 30: continue
        if paf.is_primary == False: continue 
        if paf.tname[0:5] == 'left_':
            if paf.tstart <= left_boundary_pos and paf.tend >= left_boundary_pos:
                left_anchor_paf_list.append(paf)
        elif paf.tname[0:5] == 'right':
            if paf.tstart <= right_boundary_pos and paf.tend >= right_boundary_pos:
                right_anchor_paf_list.append(paf)
        else:
            tk.eprint('ERROR! unknown template: %s' % (paf.tname))
            sys.exit(1)
    
    if len(left_anchor_paf_list) != 1 or len(right_anchor_paf_list) != 1 : return
    left_paf  = left_anchor_paf_list[0]
    right_paf  = right_anchor_paf_list[0]

    if left_paf.strand == right_paf.strand: return

    readname = left_paf.qname
    left_paf_max_repeat_size  = int((left_paf.tend - left_boundary_pos)/repeat1.repeat_unit_size) + 5
    left_paf_min_repeat_size  = tk.calculate_repeat_size_from_exact_match(left_paf.cigar, left_paf.tstart,   left_boundary_pos,  repeat1.repeat_unit_size)

    a = max(0, left_paf_min_repeat_size - 20)
    b = int(left_paf_min_repeat_size / 2.0)
    left_paf_min_repeat_size = min(a, b)
    

    initial_estimation.repeat1_count_range_dict[readname] = (left_paf_min_repeat_size, left_paf_max_repeat_size)

    
    right_paf_max_repeat_size = int((right_paf.tend - right_boundary_pos)/repeat2.repeat_unit_size) + 5
    right_paf_min_repeat_size = tk.calculate_repeat_size_from_exact_match(right_paf.cigar, right_paf.tstart, right_boundary_pos, repeat2.repeat_unit_size)
    a = max(0, right_paf_min_repeat_size - 20)
    b = int(right_paf_min_repeat_size / 2.0)
    right_paf_min_repeat_size = min(a, b)
    
    initial_estimation.repeat2_count_range_dict[readname] = (right_paf_min_repeat_size, right_paf_max_repeat_size)

    if left_paf.strand == '+': 
        candidate_start = left_paf.qstart
        candidate_end = right_paf.qlen - right_paf.qstart
        if candidate_end - candidate_start <= 0: initial_estimation.bad_reads_set.add(readname)
    else:
        candidate_start = right_paf.qstart
        candidate_end = left_paf.qlen - left_paf.qstart
        if candidate_end - candidate_start <= 0: initial_estimation.bad_reads_set.add(readname)

    initial_estimation.potential_repeat_region_dict[readname] = (candidate_start, candidate_end)
    return


def fastq_file_to_dict(in_fastq_file):

    fastq_dict = dict()

    in_fastq_f = tk.gzopen(in_fastq_file)
    while 1:
        line1 = in_fastq_f.readline()
        line2 = in_fastq_f.readline()
        line3 = in_fastq_f.readline()
        line4 = in_fastq_f.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        readname = line1.strip().split()[0][1:]
        fastq_dict[readname] = line1 + line2 + line3 + line4

    in_fastq_f.close()

    return fastq_dict

def remove_noisy_reads_2d(allele_list, ploidy, error_rate, max_mutual_overlap, max_num_components, repeat1, repeat2, in_fastq_file, out_prefix):
    tk.eprint('NOTICE: Try to remove noisy reads')
    allele_list.sort(key = lambda allele:allele.num_reads)
    num_removed_reads = 0
    while len(allele_list) > ploidy and len(allele_list) >= 2:
        if allele_list[0].num_reads * 1.5 <= allele_list[-ploidy].num_reads:
            num_removed_reads += allele_list[0].num_reads
            allele_list.pop(0)
        else:
            break

    read_repeat_joint_count_dict = dict()
    tk.eprint(f'NOTICE: There are {len(allele_list)} alleles after removing noisy reads')
    for allele in allele_list:
        for i in range(0, len(allele.readname_list)):
            readname = allele.readname_list[i]
            repeat_size1 = allele.repeat1_size_list[i]
            repeat_size2 = allele.repeat2_size_list[i]
            read_repeat_joint_count_dict[readname] = (repeat_size1, repeat_size2)
    
    remove_noisy_reads = False
    return split_alleles_using_gmm_2d (ploidy, error_rate, max_mutual_overlap, remove_noisy_reads, max_num_components, repeat1, repeat2, read_repeat_joint_count_dict, num_removed_reads, in_fastq_file, out_prefix)


def split_alleles_using_gmm_2d (ploidy, error_rate, max_mutual_overlap, remove_noisy_reads, max_num_components, repeat1, repeat2, read_repeat_joint_count_dict, num_removed_reads, in_fastq_file, out_prefix):
    
    if ploidy < 1:
        tk.eprint('ploidy must be >= 1 !\n')
        sys.exit(1)

    if len(read_repeat_joint_count_dict) < ploidy or len(read_repeat_joint_count_dict) == 1:
        tk.eprint('WARNING: No enough reads! input fastq file is: %s\n' % in_fastq_file)
        return

    cov_type = 'diag'
    dimension = 2
    probability_cutoff = 0.95

    readname_list, read_repeat_count_list = remove_outlier_reads_2d(read_repeat_joint_count_dict)
    read_repeat_count_array = np.array(read_repeat_count_list).reshape(-1, dimension)
    simulated_read_repeat_count_list = simulate_reads(read_repeat_count_list, error_rate)
    simualted_read_repeat_count_array = np.array(simulated_read_repeat_count_list).reshape(-1, dimension)

    best_n_components, final_gmm = auto_GMM_2d (simualted_read_repeat_count_array, max_num_components, cov_type, max_mutual_overlap)
    tk.eprint('NOTICE: Number of alleles = %d' % best_n_components)
    
    allele_list = create_allele_list_2d(best_n_components, final_gmm, readname_list, read_repeat_count_array, read_repeat_joint_count_dict, probability_cutoff)
    
    num_removed_reads = 0
    if remove_noisy_reads == True and len(allele_list) > ploidy:
        return remove_noisy_reads_2d(allele_list, ploidy, error_rate, max_mutual_overlap, max_num_components, repeat1, repeat2, in_fastq_file, out_prefix)
    
    allele_list.sort(key = lambda allele:allele.gmm_mean1)

    readinfo_dict = create_readinfo_dict_from_allele_list(allele_list, dimension)

    score_cut_off = calculate_log_likelyhood_cutoff(final_gmm, 0.95)

    tk.eprint('NOTICE: Writing phasing results...')
    output_phasing_results_2d(allele_list, repeat1.repeat_id, repeat2.repeat_id, in_fastq_file, out_prefix)

    tk.eprint('NOTICE: Writing to output fastq files...')
    output_phased_fastq(in_fastq_file, readinfo_dict, best_n_components, out_prefix)

    tk.eprint('NOTICE: Writing summary file...')
    output_summary_file_2d(in_fastq_file, allele_list, repeat1.repeat_id, repeat2.repeat_id, num_removed_reads, out_prefix)

    tk.eprint('NOTICE: Plotting figures...')
    plot_repeat_counts_2d(readinfo_dict, allele_list, repeat1.repeat_id, repeat2.repeat_id, out_prefix)

    scatter_plot_with_contour_2d (read_repeat_joint_count_dict, final_gmm, score_cut_off, repeat1.repeat_id, repeat2.repeat_id, allele_list, out_prefix)

    return


def calculate_log_likelyhood_cutoff(gmm, ci):
    X, labels = gmm.sample(100000)
    score_list = gmm.score_samples(X)
    X_score_list = list()
    for i in range(0, len(X)):
        X_score = [X[i], score_list[i]]
        X_score_list.append(X_score)

    X_score_list.sort(key = lambda x:x[1], reverse = True)
    cut_off_idx = int(len(X_score_list) * float(ci))
    log_likelyhood_cutoff = X_score_list[cut_off_idx][1]

    return log_likelyhood_cutoff


if __name__ == '__main__':
    main()
