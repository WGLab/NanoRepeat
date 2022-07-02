#!/usr/bin/env python3


'''
Copyright (c) 2020 Children's Hospital of Philadelphia
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
from xmlrpc.client import FastUnmarshaller, Fault
import numpy as np
import gzip
import argparse
import shutil
import random
import math
from scipy.stats import norm


import tk
from paf import *

from sklearn.mixture import GaussianMixture


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
#matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

debug = 1
        
class Readinfo:
    def __init__(self, readname):

        self.readname = readname
        self.label = -1
        self.repeat_size1 = -1
        self.repeat_size2 = -1
        self.max_prob = 0
        self.proba_list = list()
        self.qc_passed = 1

tab  = '\t' 
endl = '\n'

def parse_user_arguments():

    parser = argparse.ArgumentParser(description='Joint quantification of two adjacent tandem repeats from long-read amplicon sequencing data ')
    ### required arguments ###
    parser.add_argument('--in_fq',        required = True,  metavar = 'PATH',   type = str, help = 'input fastq file')
    parser.add_argument('--ref_fasta',    required = True,  metavar = 'PATH',   type = str, help = 'reference genome sequence in FASTA format')
    parser.add_argument('--repeat1',      required = True,  metavar = 'chr:start:end:repeat_unit:max_size', type = str, help = 'first tandem repeat, coordinates are 0-based (e.g. chr4:3074876:3074933:CAG:200)')
    parser.add_argument('--repeat2',      required = True,  metavar = 'chr:start:end:repeat_unit:max_size', type = str, help = 'second tandem repeat, coordinates are 0-based (e.g. chr4:3074946:3074966:CCG:20)')
    parser.add_argument('--out_dir',      required = True,  metavar = 'PATH',   type = str, help = 'path to the output directory')
    parser.add_argument('--version',      action='version', version='%(prog)s 0.3.0')

    ### optional arguments ### ploidy
    parser.add_argument('--num_threads',  required = False, metavar = 'INT',    type = int, default = 1,  help = 'number of threads used by minimap2 (default: 1)')
    parser.add_argument('--minimap2',     required = False, metavar = 'PATH',   type = str, default = '', help = 'path to minimap2 (default: using environment default)')
    parser.add_argument('--ploidy',       required = False, metavar = 'INT',    type = int, default = 2,  help = 'ploidy of the sample (default: 2)')
    parser.add_argument('--error_rate',   required = False, metavar = 'FLOAT',  type = float, default = 0.1,  help = 'sequencing error rate (default: 0.1)')
    parser.add_argument('--max_mutual_overlap', required = False, metavar = 'FLOAT',  type = float, default = 0.1,  help = 'max mutual overlap of two alleles in terms of repeat size distribution (default value: 0.1). If the Gaussian distribution of two alleles have more overlap than this value, the two alleles will be merged into one allele.')
    input_args = parser.parse_args()

    if input_args.ploidy < 1:
        tk.eprint('ERROR: --ploidy must be >= 1 !\n')
        sys.exit(1)

    if input_args.error_rate >= 1.0:
        tk.eprint('ERROR! --error_rate must be < 1\n')
        sys.exit(1)
    
    if input_args.max_mutual_overlap >= 1.0:
        tk.eprint('ERROR! --max_mutual_overlap must be < 1\n')
        sys.exit(1)

    tk.check_input_file_exists(input_args.in_fq)
    tk.check_input_file_exists(input_args.ref_fasta)

    input_args.in_fq       = os.path.abspath(input_args.in_fq)
    input_args.ref_fasta   = os.path.abspath(input_args.ref_fasta)
    input_args.out_dir     = os.path.abspath(input_args.out_dir)

    if input_args.minimap2 == '': 
        input_args.minimap2 = tk.find_executable_path('minimap2')
        if not input_args.minimap2:
            tk.eprint('ERROR! minimap2 was not found! Please supply the path to --minimap2')
            sys.exit(1)
        else:
            if os.path.exists(input_args.minimap2):
                tk.eprint('NOTICE: found path to minimap2: %s' % input_args.minimap2)
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

    tk.create_dir(input_args.out_dir)

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

    in_fastq_prefix = os.path.splitext(os.path.split(input_args.in_fq)[1])[0]
    temp_out_dir = os.path.join(input_args.out_dir, '%s.NanoRepeat.temp' % in_fastq_prefix)
    tk.create_dir(temp_out_dir)

    repeat_chrom_seq = tk.read_one_chr_from_fasta_file(input_args.ref_fasta, repeat1.chrom)
    if len(repeat_chrom_seq) == 0:
        tk.eprint('ERROR: ref_fasta file: %s has no valid sequence!\n' % input_args.ref_fasta)
        sys.exit(1)
    
    platform = 'ont'

    tk.eprint(f'NOTICE: input file is: {input_args.in_fq }')
    initial_estimation = initial_estimate_repeat_size(input_args.minimap2, repeat_chrom_seq, input_args.in_fq, platform, input_args.num_threads, repeat1, repeat2, max_anchor_len, temp_out_dir)

    final_estimation = fine_tune_read_count(initial_estimation, input_args.in_fq, repeat_chrom_seq, repeat1, repeat2, input_args.minimap2, platform, input_args.num_threads, temp_out_dir)
    
    jointly_split_alleles_using_gmm (input_args.ploidy, input_args.error_rate, input_args.max_mutual_overlap, repeat1, repeat2, final_estimation, input_args.in_fq, input_args.out_dir)

    tk.eprint('NOTICE: program finished. Output files are here: %s\n' % input_args.out_dir)

    shutil.rmtree(temp_out_dir)

    return


def output_repeat_size_file(read_repeat_joint_count_dict, readinfo_dict, in_fastq_file, repeat_region_list, repeat_size_file, ploidy):

    repeat_size_fp = open(repeat_size_file, 'w')
    repeat_size_fp.write('##input_fastq=%s\n' % in_fastq_file)
    repeat_size_fp.write('#readname\t%s\t%s\tallele_id\tprobability\tQC\n' % (repeat_region_list[0].repeat_id, repeat_region_list[1].repeat_id))
    for readname in read_repeat_joint_count_dict:
        repeat_size1, repeat_size2 = read_repeat_joint_count_dict[readname]
        if readname in readinfo_dict:
            readinfo = readinfo_dict[readname]
            allele_id = readinfo.label + 1
            if readinfo.qc_passed:
                qc = 'PASS'
            else:
                qc = 'FAIL'
            max_prob = readinfo.max_prob
            repeat_size_fp.write('%s\t%d\t%d\t%d\t%.8f\t%s' % (readname, repeat_size1, repeat_size2, allele_id, max_prob, qc))
            for prob in readinfo.proba_list:
                repeat_size_fp.write('\t%.8f' % prob )
            repeat_size_fp.write('\n')

        else:
            allele_id = 'N.A.'
            max_prob = 'N.A.'
            qc = 'FAIL'
            repeat_size_fp.write('%s\t%d\t%d\t%s\t%s\t%s' % (readname, repeat_size1, repeat_size2, allele_id, max_prob, qc))
            for i in range(0, ploidy):
                repeat_size_fp.write('\tN.A.')
            repeat_size_fp.write('\n')
    repeat_size_fp.close()

    return

def fine_tune_read_count(initial_estimation, in_fastq_file, repeat_chrom_seq, repeat1, repeat2, minimap2, platform, num_threads, out_dir):

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

    round2_estimation = round2_estimation_of_repeat_size(initial_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, platform, num_threads, out_dir)
    
    if round2_estimation.step_size1 > 1 and round2_estimation.step_size2 > 1:
        final_estimation = round3_estimation_of_repeat_size(initial_estimation, round2_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, platform, num_threads, out_dir)
    else:
        final_estimation = round2_estimation
    
    return final_estimation

def round3_estimation_of_repeat_size(initial_estimation, round2_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, platform, num_threads, out_dir):

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

    preset = tk.get_preset_for_minimap2(platform)
    round3_paf_file =  os.path.join(out_dir, 'round3.paf')
    round3_paf_fp = open(round3_paf_file, 'w')
    round3_paf_fp.close()

    for repeat_count1 in range(min_size1, max_size1):
        for repeat_count2 in range(min_size2, max_size2):
            tmp_fastq_file = os.path.join(out_dir, '%d.%d.round3.fastq' % (repeat_count1, repeat_count2) )
            tmp_fastq_fp = open(tmp_fastq_file, 'w')
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

                tmp_fastq_fp.write(fastq_dict[readname])
                n_tmp_reads += 1
            tmp_fastq_fp.close()
            if n_tmp_reads == 0: 
                tk.rm(tmp_fastq_file)
                continue
            tmp_ref_file = os.path.join(out_dir, '%d.%d.ref.fasta' % (repeat_count1, repeat_count2))
            build_fasta_template_for_two_repeats(left_anchor_seq, mid_anchor_seq, right_anchor_seq, repeat1, repeat2, repeat_count1, repeat_count2, tmp_ref_file)
            cmd = f'{minimap2} -A 2 -B 6 -c --eqx -t {num_threads} -x {preset} {tmp_ref_file} {tmp_fastq_file} >> {round3_paf_file} 2> /dev/null'
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

def round2_estimation_of_repeat_size(initial_estimation, fastq_dict, repeat_chrom_seq, repeat1, repeat2, minimap2, platform, num_threads, out_dir):
    
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

    preset = tk.get_preset_for_minimap2(platform)
    round2_paf_file =  os.path.join(out_dir, 'round2.paf')
    round2_paf_fp = open(round2_paf_file, 'w')
    round2_paf_fp.close()

    for repeat_count1 in range(repeat1.round1_min_size, repeat1.round1_max_size + 1, step_size1):
        for repeat_count2 in range(repeat2.round1_min_size, repeat2.round1_max_size + 1, step_size2):
            tmp_fastq_file = os.path.join(out_dir, '%d.%d.round2.fastq' % (repeat_count1, repeat_count2) )
            tmp_fastq_fp = open(tmp_fastq_file, 'w')
            n_tmp_reads = 0
            for readname in fastq_dict:
                if readname not in initial_estimation.repeat1_count_range_dict: continue
                if readname not in initial_estimation.repeat2_count_range_dict: continue
                min1, max1 = initial_estimation.repeat1_count_range_dict[readname]
                min2, max2 = initial_estimation.repeat2_count_range_dict[readname]
                if repeat_count1 >= min1 and repeat_count1 < max1 and repeat_count2 >= min2 and repeat_count2 < max2:                
                    tmp_fastq_fp.write(fastq_dict[readname])
                    n_tmp_reads += 1
            tmp_fastq_fp.close()
            if n_tmp_reads == 0:
                tk.rm(tmp_fastq_file)
                continue
            tmp_ref_file = os.path.join(out_dir, '%d.%d.ref.fasta' % (repeat_count1, repeat_count2))
            build_fasta_template_for_two_repeats(left_anchor_seq, mid_anchor_seq, right_anchor_seq, repeat1, repeat2, repeat_count1, repeat_count2, tmp_ref_file)
            cmd = f'{minimap2} -A 2 -B 6 -c --eqx -t {num_threads} -x {preset} {tmp_ref_file} {tmp_fastq_file} >> {round2_paf_file} 2> /dev/null'
            tk.run_system_cmd(cmd)
            tk.rm(tmp_ref_file)
            tk.rm(tmp_fastq_file)

    round2_estimation = estimate_two_repeats_from_paf(round2_paf_file, left_anchor_len, mid_anchor_len, repeat1, repeat2)
    round2_estimation.step_size1 = step_size1
    round2_estimation.step_size2 = step_size2

    return round2_estimation

def estimate_two_repeats_from_paf(in_paf_file, left_anchor_len, mid_anchor_len, repeat1, repeat2):

    repeat_size_estimation = RepeatSize()

    in_paf_fp = open(in_paf_file, 'r')
    read_align_score_dict = dict()
    while 1:
        line = in_paf_fp.readline()
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

    in_paf_fp.close()

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

    tmp_ref_fp = open(tmp_ref_file, 'w')
    tmp_ref_fp.write('>%d-%d\n' % (repeat_count1, repeat_count2))
    seq = left_anchor_seq + repeat1.repeat_unit * repeat_count1 + mid_anchor_seq + repeat2.repeat_unit * repeat_count2 + right_anchor_seq + '\n'
    tmp_ref_fp.write(seq)
    tmp_ref_fp.close()

    return

def initial_estimate_repeat_size(minimap2, repeat_chrom_seq, in_fastq_file, platform, num_threads, repeat1, repeat2, max_anchor_len, out_dir):

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

    left_template_fasta_fp = open(left_template_fasta_file, 'w')
    left_template_fasta_fp.write('>%s\n' % left_template_name)
    left_template_fasta_fp.write('%s\n' % left_template_seq)
    left_template_fasta_fp.close()


    right_template_fasta_file = os.path.join(out_dir, 'round1_templates.right_anchor.fasta')
    right_template_name = 'right_anchor_%d_%d_%s_revc' % (len(right_anchor_seq), repeat2.max_size, repeat2.repeat_unit)
    right_template_seq  = repeat2.repeat_unit * repeat2.max_size + right_anchor_seq
    right_template_seq  = tk.rev_comp(right_template_seq)
    
    right_template_fasta_fp = open(right_template_fasta_file, 'w')
    right_template_fasta_fp.write('>%s\n' % right_template_name)
    right_template_fasta_fp.write('%s\n' % right_template_seq)
    right_template_fasta_fp.close()

    round1_paf_file = os.path.join(out_dir, 'round1.paf')
    round1_left_anchor_paf_file  = os.path.join(out_dir, 'round1.left.paf')
    round1_right_anchor_paf_file = os.path.join(out_dir, 'round1.right.paf')
    preset = tk.get_preset_for_minimap2(platform)

    cmd = f'{minimap2} -A 2 -B 6 -c --eqx -t {num_threads} -x {preset} {left_template_fasta_file} {in_fastq_file} > {round1_left_anchor_paf_file} 2> /dev/null'
    tk.run_system_cmd(cmd)
    
    cmd = f'{minimap2} -A 2 -B 6 -c --eqx -t {num_threads} -x {preset} {right_template_fasta_file} {in_fastq_file} > {round1_right_anchor_paf_file} 2> /dev/null'
    tk.run_system_cmd(cmd)

    cmd = f'cat {round1_left_anchor_paf_file} {round1_right_anchor_paf_file} | sort -k1 - > {round1_paf_file}'
    tk.run_system_cmd(cmd)

    tk.rm(round1_left_anchor_paf_file)
    tk.rm(round1_right_anchor_paf_file)
    initial_estimation = round1_estimation_from_paf(round1_paf_file, repeat1, repeat2, left_anchor_len, right_anchor_len)

    return initial_estimation

def round1_estimation_from_paf(round1_paf_file, repeat1, repeat2, left_anchor_len, right_anchor_len):

    initial_estimation = Round1Estimation()
    round1_paf_fp = open(round1_paf_file, 'r')
    read_paf_list = list()
    while 1:
        line = round1_paf_fp.readline()
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

    round1_paf_fp.close()
    
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

    in_fastq_fp = tk.gzopen(in_fastq_file)
    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        readname = line1.strip().split()[0][1:]
        fastq_dict[readname] = line1 + line2 + line3 + line4

    in_fastq_fp.close()

    return fastq_dict


def jointly_split_alleles_using_gmm (ploidy, error_rate, max_mutual_overlap, repeat1, repeat2, final_estimation, in_fastq_file, out_dir):
    
    in_fastq_prefix  = os.path.splitext(os.path.split(in_fastq_file)[1])[0]
    out_prefix = os.path.join(out_dir, '%s.JointGMM' % (in_fastq_prefix))
    out_fastq_prefix = os.path.join(out_dir, '%s.JointGMM.qc_passed' % (in_fastq_prefix))
    repeat_size_file = os.path.join(out_dir, '%s.repeat_size.txt' % in_fastq_prefix)

    read_repeat_joint_count_dict = dict()
    for readname in final_estimation.repeat1_count_dict:
        if readname not in final_estimation.repeat2_count_dict: continue
        size1 = final_estimation.repeat1_count_dict[readname]
        size2 = final_estimation.repeat2_count_dict[readname]
        read_repeat_joint_count_dict[readname] = (size1, size2)

    repeat_region_list = [repeat1, repeat2]

    if len(read_repeat_joint_count_dict) < ploidy or len(read_repeat_joint_count_dict) == 1:
        tk.eprint('WARNING: No enough reads! input fastq file is: %s\n' % in_fastq_file)
        readinfo_dict = dict()
        output_repeat_size_file(read_repeat_joint_count_dict, readinfo_dict, in_fastq_file, repeat_region_list, repeat_size_file, ploidy)
        return

    if ploidy < 1:
        tk.eprint('ploidy must be >= 1 !\n')
        sys.exit(1)

    proba_cutoff = 0.95
    #cov_type = 'tied'
    cov_type = 'diag'
    dimension = 2
    min_count1, max_count1, min_count2, max_count2  = analysis_outlier_2d(read_repeat_joint_count_dict)

    readname_list = list() # readnames after removing outliers
    read_repeat_count_list = list()

    for readname in read_repeat_joint_count_dict:
        repeat_count1, repeat_count2 = read_repeat_joint_count_dict[readname]
        if repeat_count1 < min_count1 or repeat_count1 > max_count1: continue 
        if repeat_count2 < min_count2 or repeat_count2 > max_count2: continue 
        readname_list.append(readname)
        read_repeat_count_list.append(repeat_count1)
        read_repeat_count_list.append(repeat_count2)

    read_repeat_count_array = np.array(read_repeat_count_list).reshape(-1, dimension)

    min_std = 1.0
    simulated_read_repeat_count_list = read_repeat_count_list * 50
    for i in range(0, len(simulated_read_repeat_count_list)):
        std = min_std + 1.25 * error_rate/3.0 * simulated_read_repeat_count_list[i]
        random_error = random.gauss(0, std)
        simulated_read_repeat_count_list[i] += random_error


    simualted_read_repeat_count_array = np.array(simulated_read_repeat_count_list).reshape(-1, dimension)

    best_n_components, final_gmm = auot_GMM (simualted_read_repeat_count_array, ploidy+2, cov_type, max_mutual_overlap)
    tk.eprint('NOTICE: number of alleles = %d' % best_n_components)
    old_read_label_list = list(final_gmm.predict(read_repeat_count_array))
    proba2darray = final_gmm.predict_proba(read_repeat_count_array)

    repeat_region_list[0].is_main = 1
    main_repeat_idx = 0
    for i in range(0, len(repeat_region_list)):
        if repeat_region_list[i].is_main:
            main_repeat_idx = i
            break

    old_cluster_mean_list = list()
    for means in final_gmm.means_:
        old_cluster_mean_list.append(means[main_repeat_idx])

    read_label_list, old_label_to_new_label_dict, new_label_to_old_label_dict = sort_label_by_cluster_mean(old_read_label_list, old_cluster_mean_list)

    readinfo_dict = dict()
    for i in range(0, len(readname_list)):
        readname = readname_list[i]
        repeat_count1, repeat_count2 = read_repeat_joint_count_dict[readname]
        read_label = read_label_list[i]
        max_prob = max(proba2darray[i])
        readinfo = Readinfo(readname)
        readinfo.label = read_label
        readinfo.max_prob = max_prob
        readinfo.proba_list = list(proba2darray[i])
        readinfo.repeat_size1 = repeat_count1
        readinfo.repeat_size2 = repeat_count2
        readinfo_dict[readname] = readinfo
    
    each_allele_repeat_count_3d_list = [0] * best_n_components # each_allele_repeat_count_3d_list[label][repeat_id] = list of repeat size
    for label in range(0, len(each_allele_repeat_count_3d_list)):
        each_allele_repeat_count_3d_list[label] = [0] * 2
        for repeat_id in range(0, 2):
            each_allele_repeat_count_3d_list[label][repeat_id] = list()

    for readname in readinfo_dict:
        readinfo = readinfo_dict[readname]
        each_allele_repeat_count_3d_list[readinfo.label][0].append(readinfo.repeat_size1)
        each_allele_repeat_count_3d_list[readinfo.label][1].append(readinfo.repeat_size2)

    allele_predicted_repeat_size_2d_list = [0] * best_n_components
    for label in range(0, len(allele_predicted_repeat_size_2d_list)):
        allele_predicted_repeat_size_2d_list[label] = [0] * 2
        allele_predicted_repeat_size_2d_list[label][0] = int(np.median(each_allele_repeat_count_3d_list[label][0]) + 0.5)
        allele_predicted_repeat_size_2d_list[label][1] = int(np.median(each_allele_repeat_count_3d_list[label][1]) + 0.5)

    score_cut_off = calculate_log_likelyhood_cutoff(final_gmm, 0.99)
    label_qc_failed_reads(readinfo_dict, final_gmm, proba_cutoff, score_cut_off)

    tk.eprint('NOTICE: writing to repeat size file.')
    output_repeat_size_file(read_repeat_joint_count_dict, readinfo_dict, in_fastq_file, repeat_region_list, repeat_size_file, ploidy)
    tk.eprint('NOTICE: writing to output fastq files.')
    joint_gmm_output_fastq(in_fastq_file, readinfo_dict, best_n_components, out_fastq_prefix)

    tk.eprint('NOTICE: writing to output summary file.')
    joint_gmm_output_summary_file(in_fastq_file, readinfo_dict, allele_predicted_repeat_size_2d_list, repeat_region_list, out_prefix)

    tk.eprint('NOTICE: plotting figures.')
    joint_gmm_plot_repeat_counts(readinfo_dict, allele_predicted_repeat_size_2d_list, repeat_region_list, out_prefix)

    joint_gmm_scatter_plot_with_contour (read_repeat_joint_count_dict, final_gmm, score_cut_off, repeat_region_list, out_prefix)

    return

def joint_gmm_scatter_plot_with_contour (read_repeat_joint_count_dict, final_gmm, score_cut_off, repeat_region_list, out_prefix):

    scatter_plot_file = out_prefix + '.scatter.png'

    xlabel = '%s repeat size' % (repeat_region_list[0].repeat_id)
    ylabel = '%s repeat size' % (repeat_region_list[1].repeat_id)

    X = list()
    Y = list()
    
    for readname in read_repeat_joint_count_dict:
        x, y = read_repeat_joint_count_dict[readname]
        X.append(x)
        Y.append(y)

    X = np.array(X)
    Y = np.array(Y)

    X, Y, Z= countxy(X, Y)
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 4)

    try:
        ax.scatter(X, Y, c=Z, s=15, edgecolor='')
    except:
        ax.scatter(X, Y, c=Z, s=15, edgecolor=['none'])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    norm = Normalize(vmin = np.min(Z), vmax = np.max(Z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Count')
    
    xmin = min(X)
    xmax = max(X)
    ymin = min(Y)
    ymax = max(Y)

    figure_max_x = int(((xmax * 1.4)/10.0 +1 )* 10)
    figure_max_y = int(((ymax * 1.4)/10.0 +1 )* 10)

    figure_min_x = int(((xmin / 1.4)/10.0 -1 )* 10)
    figure_min_y = int(((ymin / 1.4)/10.0 -1 )* 10)

    if figure_min_x < 10: figure_min_x = 0
    if figure_min_y < 10: figure_min_y = 0

    plt.xlim(figure_min_x, figure_max_x)
    plt.ylim(figure_min_y, figure_max_y)

    a = np.linspace(figure_min_x, figure_max_x, 400)
    b = np.linspace(figure_min_y, figure_max_y, 400)

    A, B = np.meshgrid(a, b)
    AA = np.array([A.ravel(), B.ravel()]).T
    C = final_gmm.score_samples(AA)
    C = C.reshape(A.shape)

    CS = plt.contour(A, B, C, levels=[score_cut_off], linestyles = 'dashed', colors = 'grey')

    plt.savefig(scatter_plot_file, dpi = 300)
    plt.close('all')
    return

def countxy(x, y):
    count_dict = dict()
    for i in range(0, len(x)):
        key = '%d\t%d' % (x[i], y[i])
        if key not in count_dict:
            count_dict[key] = 1
        else:
            count_dict[key] += 1
    x_list = list()
    y_list = list()
    z_list = list()
    for key in count_dict:
        x, y = key.split('\t')
        x = int(x)
        y = int(y)
        z = count_dict[key]
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)
    
    return x_list, y_list, z_list

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

def joint_gmm_plot_repeat_counts(readinfo_dict, allele_predicted_repeat_size_2d_list, repeat_region_list, out_prefix):

    
    hist2d_figure_file = out_prefix + '.hist2d.png'
    repeat1_hist_figure_file  = out_prefix + '.%s.hist.png' % (repeat_region_list[0].repeat_id)
    repeat2_hist_figure_file  = out_prefix + '.%s.hist.png' % (repeat_region_list[1].repeat_id)

    num_alleles = len(allele_predicted_repeat_size_2d_list)
    x_list = list()
    y_list = list()

    x_2d_list = [0] * num_alleles
    y_2d_list = [0] * num_alleles
    for i in range(0, num_alleles):
        x_2d_list[i] = list()
        y_2d_list[i] = list()
    for readname in readinfo_dict:
        readinfo = readinfo_dict[readname]
        if readinfo.qc_passed == 0: continue

        x = readinfo.repeat_size1
        y = readinfo.repeat_size2
        x_list.append(x)
        y_list.append(y)

        x_2d_list[readinfo.label].append(x)
        y_2d_list[readinfo.label].append(y)

    xmin = int(min(x_list))
    xmax = int(max(x_list))+1
    ymin = int(min(y_list))
    ymax = int(max(y_list))+1

    if xmax - xmin <= 200:
        b1 = range(xmin - 1, xmax + 2)
    else:
        b1 = range(xmin - 1, xmax + 2, int(float(xmax - xmin)/200.0 + 0.5))

    if ymax - ymin <= 200:
        b2 = range(ymin - 1, ymax + 2)
    else:
        b2 = range(ymin - 1, ymax + 2, int(float(ymax - ymin)/200.0 + 0.5))
    
    repeat1_predicted_size_list = list()
    repeat2_predicted_size_list = list()
    for label in range(0, len(allele_predicted_repeat_size_2d_list)):
        repeat1_predicted_size_list.append(allele_predicted_repeat_size_2d_list[label][0])
        repeat2_predicted_size_list.append(allele_predicted_repeat_size_2d_list[label][1])

    plot_hist2d(x_list, y_list, b1, b2, repeat_region_list[0].repeat_id, repeat_region_list[1].repeat_id, hist2d_figure_file)
    plot_hist1d(x_2d_list, b1, repeat_region_list[0].repeat_id, repeat1_predicted_size_list, repeat1_hist_figure_file)
    plot_hist1d(y_2d_list, b2, repeat_region_list[1].repeat_id, repeat2_predicted_size_list, repeat2_hist_figure_file)

    return 

def plot_hist1d(x_2d_list, b, repeat_id, predicted_size_list, out_file):

    plt.figure (figsize=(6, 4))
    
    max_x = 0
    for x in x_2d_list:
        if len(x) == 0: continue
        plt.hist(x, bins = b)
        if max_x < max(x):
            max_x = max(x)

    for repeat_size in predicted_size_list:
        plt.axvline(x=repeat_size+0.5, color = 'grey', linestyle = ':')


    plt.title('Repeat size distribution (%s)' % repeat_id)
    plt.xlabel('repeat size')
    plt.ylabel('number of reads')
    
    max_x = int(((max_x * 1.618)/5.0 +1 )* 5)
    plt.xlim(0, max_x)

    plt.savefig(out_file, dpi=300)
    plt.close('all')

    return

def plot_hist2d(x_list, y_list, b1, b2, repeat_id1, repeat_id2, out_file):

    plt.figure (figsize=(8, 4))
    plt.hist2d (x_list, y_list, [b1, b2], cmap = 'binary')
    plt.colorbar()
    plt.title('2D histogram of repeat size')
    plt.xlabel('repeat size (%s)' % repeat_id1)
    plt.ylabel('repeat size (%s)' % repeat_id2)
    xmin = min(x_list)
    xmax = max(x_list)

    if debug:
        plt.xlim(xmin, xmax)
        plt.ylim(0, (xmax-xmin)/2)

    plt.savefig(out_file, dpi=300)
    plt.close('all')

    return

def joint_gmm_output_summary_file(in_fastq_file, readinfo_dict, allele_predicted_repeat_size_2d_list, repeat_region_list, out_prefix):

    out_summray_file = out_prefix + '.summary.txt'
    out_summray_fp = open(out_summray_file, 'w')
    summary_header = '#input_fastq\tmethod'
    summary_info = '%s\tJointGMM' % (in_fastq_file)

    # allele_predicted_repeat_size_2d_list[label][repeat_id] = predicted_repeat_size
    num_alleles = len(allele_predicted_repeat_size_2d_list)
    summary_header += '\tnum_alleles'
    summary_info   += '\t%d' % (num_alleles)

    allele_num_reads_list = [0] * num_alleles
    for readname in readinfo_dict:
        readinfo = readinfo_dict[readname]
        if readinfo.qc_passed == 0: continue
        allele_num_reads_list[readinfo.label] += 1

    for label in range(0, num_alleles):
        allele_id = label + 1
        summary_header += '\tallele%d_num_reads' % (allele_id)
        summary_info   += '\t%d' % allele_num_reads_list[label]
        for i in range(0, 2):
            summary_header += '\t%s_repeat_size%d' % (repeat_region_list[i].repeat_id, allele_id)
            summary_info   += '\t%d' % (allele_predicted_repeat_size_2d_list[label][i])

    out_summray_fp.write(summary_header + '\n')
    out_summray_fp.write(summary_info + '\n')

    return
    
def joint_gmm_output_fastq(in_fastq_file, readinfo_dict, num_alleles, out_prefix):

    out_allele_fastq_file_list = list()
    for label in range(0, num_alleles):
        allele_id = label + 1
        out_allele_fastq_file = out_prefix + '.allele%d.fastq' % (allele_id)
        out_allele_fastq_file_list.append(out_allele_fastq_file)

    out_allele_fastq_fp_list = list()
    for i in range(0, len(out_allele_fastq_file_list)):
        out_allele_fastq_fp = open(out_allele_fastq_file_list[i], 'w')
        out_allele_fastq_fp_list.append(out_allele_fastq_fp)
    
    if '.gz' == in_fastq_file[-3:]:
        in_fastq_fp = gzip.open(in_fastq_file, 'rt')
    else:
        in_fastq_fp = open(in_fastq_file, 'rt')

    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        readname = line1.strip().split()[0][1:]
        if readname not in readinfo_dict: continue
        if readinfo_dict[readname].qc_passed == 0: continue
    
        label = readinfo_dict[readname].label
        out_allele_fastq_fp_list[label].write(line1 + line2 + line3 + line4)

    in_fastq_fp.close()

    for i in range(0, len(out_allele_fastq_fp_list)):
        out_allele_fastq_fp_list[i].close()

    return


def label_qc_failed_reads(readinfo_dict, final_gmm, proba_cutoff, score_cut_off):

    for readname in readinfo_dict:
        readinfo_dict[readname].qc_passed = 1
        continue
        readinfo = readinfo_dict[readname]
        qc_failed = 0

        x = [readinfo.repeat_size1, readinfo.repeat_size2]
        x = np.array(x).reshape(1, 2)
        score = final_gmm.score_samples(x)
        if score[0] < score_cut_off: qc_failed = 1

        if readinfo.max_prob < proba_cutoff: qc_failed = 1

        if qc_failed: readinfo_dict[readname].qc_passed = 0

    return


def auot_GMM(X, max_num_components, cov_type, max_mutual_overlap):
    for n in range(max_num_components, 0, -1):
        gmm = GaussianMixture(n_components=n, covariance_type=cov_type, n_init=10).fit(X)
        if n == 1: 
            return n, gmm
        too_close = False
        for i in range(0, n):
            for j in range(i+1, n):
                component_i_mean1 = gmm.means_[i][0]
                component_i_mean2 = gmm.means_[i][1]
                component_j_mean1 = gmm.means_[j][0]
                component_j_mean2 = gmm.means_[j][1]
   
                component_i_cov1 = gmm.covariances_[i][0]
                component_i_cov2 = gmm.covariances_[i][1]
                component_j_cov1 = gmm.covariances_[j][0]
                component_j_cov2 = gmm.covariances_[j][1]
                
                component_i_sd1 = max(1.0, math.sqrt(component_i_cov1))
                component_i_sd2 = max(1.0, math.sqrt(component_i_cov2))
                component_j_sd1 = max(1.0, math.sqrt(component_j_cov1))
                component_j_sd2 = max(1.0, math.sqrt(component_j_cov2))

                a1 = 1.0-max_mutual_overlap
                a2 = max_mutual_overlap
                interval_i1 = (norm.isf(a1, component_i_mean1, component_i_sd1), norm.isf(a2, component_i_mean1, component_i_sd1))
                interval_i2 = (norm.isf(a1, component_i_mean2, component_i_sd2), norm.isf(a2, component_i_mean2, component_i_sd2))
                interval_j1 = (norm.isf(a1, component_j_mean1, component_j_sd1), norm.isf(a2, component_j_mean1, component_j_sd1))
                interval_j2 = (norm.isf(a1, component_j_mean2, component_j_sd2), norm.isf(a2, component_j_mean2, component_j_sd2))
                if has_overlap(interval_i1, interval_j1) and has_overlap(interval_i2, interval_j2):
                    too_close = True
                    break

        if not too_close:
            return n, gmm
    

def has_overlap(interval1, interval2):
    start = min(interval1[1], interval2[1])
    end = max(interval1[0], interval2[0])
    if end - start <= 0: 
        return True
    else:
        return False

def analysis_outlier(read_repeat_count_dict):

    read_repeat_count_list = list()
    for readname in read_repeat_count_dict:
        repeat_count = read_repeat_count_dict[readname]
        read_repeat_count_list.append(repeat_count)

    min_repeat_count_cutoff, max_repeat_count_cutoff = get_outlier_cutoff_from_list(read_repeat_count_list)
 
    return min_repeat_count_cutoff, max_repeat_count_cutoff

def analysis_outlier_2d(read_repeat_joint_count_dict):

    repeat_count1_list = list()
    repeat_count2_list = list()

    for readname in read_repeat_joint_count_dict:
        repeat_count1, repeat_count2 = read_repeat_joint_count_dict[readname]
        repeat_count1_list.append(repeat_count1)
        repeat_count2_list.append(repeat_count2)

    min_count1, max_count1 = get_outlier_cutoff_from_list(repeat_count1_list)
    min_count2, max_count2 = get_outlier_cutoff_from_list(repeat_count2_list)
    
    return min_count1, max_count1, min_count2, max_count2

def get_outlier_cutoff_from_list(repeat_count_list):

    mean = np.mean(repeat_count_list)
    std = np.std(repeat_count_list)

    min_repeat_count_cutoff = mean - 3 * std
    if min_repeat_count_cutoff < 0: min_repeat_count_cutoff = 0
    max_repeat_count_cutoff = mean + 3 * std

    return min_repeat_count_cutoff, max_repeat_count_cutoff

def sort_label_by_cluster_mean(old_read_label_list, cluster_mean_list):

    l = list()
    for i in range(0, len(cluster_mean_list)):
        label = i
        cluster_mean = cluster_mean_list[i]
        l.append((label, cluster_mean))

    l = sorted(l, key=lambda x:x[1])

    old_label_to_new_label_dict = dict()
    new_label_to_old_label_dict = dict()

    for i in range(0, len(l)):
        new_label = i
        old_label = l[i][0]
        old_label_to_new_label_dict[old_label] = new_label
        new_label_to_old_label_dict[new_label] = old_label

    new_read_label_list = list()

    for i in range(0, len(old_read_label_list)):
        old_label = old_read_label_list[i]
        new_label = old_label_to_new_label_dict[old_label]
        new_read_label_list.append(new_label)

    return new_read_label_list, old_label_to_new_label_dict, new_label_to_old_label_dict


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
        self.repeat_id = '_'.join(col_list[0:4])
        #self.start -= 1
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


if __name__ == '__main__':
    main()
