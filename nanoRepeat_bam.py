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
import string
import sys
import shutil
import numpy as np
import tk
from repeat_region import *
from paf import *


from numpy.f2py.auxfuncs import containsmodule
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import List

def extract_ref_sequence(ref_fasta_dict, repeat_region:RepeatRegion):

    MIN_ANCHOR_LEN = 10
    if repeat_region.anchor_len < MIN_ANCHOR_LEN:
        tk.eprint('WARNING: anchor_len must be >= 10!')
        repeat_region.anchor_len = 10

    chr_name = repeat_region.chrom
    if chr_name not in ref_fasta_dict:
        if 'chr' == chr_name[0:3]:
            chr_name = chr_name[3:]
        else:
            chr_name = 'chr' + chr_name
    if chr_name not in ref_fasta_dict:
        tk.eprint(f'ERROR! the following chromosome in repeat region bed file was not found in the reference fasta file: {repeat_region.chrom}')
        sys.exit()
    
    chr_seq = ref_fasta_dict[chr_name]
    chr_len = len(chr_seq)

    if repeat_region.start_pos > chr_len:
        tk.eprint(f'ERROR! the repeat start position is larger than chromosome length: {repeat_region.start_pos}. chromosome length = {chr_len}')
        sys.exit()
    
    if repeat_region.start_pos < 0:
        tk.eprint(f'ERROR! the repeat start position < 0')
        sys.exit()

    if repeat_region.end_pos > chr_len + 1:
        tk.eprint(f'ERROR! the repeat end position is larger than chromosome length: {repeat_region.end_pos}. chromosome length = {chr_len}')
        sys.exit()
    
    if repeat_region.end_pos < repeat_region.start_pos:
        tk.eprint(f'ERROR! end position is smaller than start position: {repeat_region.to_invertal()}')
        sys.exit()

    start_pos = repeat_region.start_pos - repeat_region.anchor_len
    end_pos = repeat_region.end_pos + repeat_region.anchor_len
    if start_pos < 0: 
        start_pos = 0
    
    if end_pos > chr_len:
        end_pos = chr_len
    
    repeat_region.left_anchor_seq  = chr_seq[start_pos:repeat_region.start_pos]
    repeat_region.left_anchor_len  = len(repeat_region.left_anchor_seq)

    repeat_region.right_anchor_seq = chr_seq[repeat_region.end_pos:end_pos]
    repeat_region.right_anchor_len = len(repeat_region.right_anchor_seq)

    repeat_region.mid_ref_seq      = chr_seq[repeat_region.start_pos:repeat_region.end_pos] # repeat sequence in the ref genome

    if repeat_region.left_anchor_len == 0 and repeat_region.right_anchor_len == 0:
        tk.eprint('ERROR! there is no flanking sequence around the repeat region! ')
        sys.exit(1)

    if repeat_region.left_anchor_len < MIN_ANCHOR_LEN and repeat_region.right_anchor_len < MIN_ANCHOR_LEN:
        tk.eprint(f'ERROR! both left and right flanking sequences are less than {MIN_ANCHOR_LEN} bp. The quantification might not be accurate!')
        sys.exit(1)
    
    return

def check_anchor_mapping(one_read_anchor_paf_list: List[PAF]):

    if len(one_read_anchor_paf_list) == 0:
        return False
    
    if len(one_read_anchor_paf_list) == 1:
        return True
    
    if one_read_anchor_paf_list[0].align_len < 10: 
        return False
    
    if one_read_anchor_paf_list[0].align_score > 1.5 * one_read_anchor_paf_list[1].align_score and one_read_anchor_paf_list[0].mapq > 30:
        return True
    else:
        return False

def find_anchor_locations_for1read(read_paf_list: List[PAF], repeat_region: RepeatRegion):

    left_anchor_paf_list  = []
    right_anchor_paf_list = []
    for paf in read_paf_list:
        if paf.tname.startswith('left_anchor'):
            left_anchor_paf_list.append(paf)
        elif paf.tname.startswith('right_anchor'):
            right_anchor_paf_list.append(paf)

    left_anchor_paf_list.sort(key = lambda paf:paf.align_score, reverse=True)
    right_anchor_paf_list.sort(key = lambda paf:paf.align_score, reverse=True)
    read = Read()
    read.read_name = read_paf_list[0].qname
    read.full_read_len = read_paf_list[0].qlen
    read.both_anchors_are_good = False
    read.left_anchor_is_good = check_anchor_mapping(left_anchor_paf_list)
    read.right_anchor_is_good = check_anchor_mapping(right_anchor_paf_list)
    
    if read.left_anchor_is_good == False: return
    if read.right_anchor_is_good == False: return
    
    repeat_region_length = 0
    read.left_anchor_paf = left_anchor_paf_list[0]
    read.right_anchor_paf = right_anchor_paf_list[0]

    # Important: if strand == '-', all positions are based on the reverse strand
    if read.left_anchor_is_good and read.right_anchor_is_good and read.left_anchor_paf.strand == read.right_anchor_paf.strand:
        repeat_region_length = read.right_anchor_paf.qstart - read.left_anchor_paf.qend

    if repeat_region_length > -10:
        read.both_anchors_are_good = True
        read.dist_between_anchors = repeat_region_length
    
    if read.both_anchors_are_good == False: return

    repeat_region.read_dict[read.read_name] = read

    repeat_region.buffer_len = 100
    if read.both_anchors_are_good:
        read.core_seq_start_pos = read.left_anchor_paf.qend - repeat_region.buffer_len
        read.core_seq_end_pos   = read.right_anchor_paf.qstart + repeat_region.buffer_len
        read.mid_seq_start_pos  = read.left_anchor_paf.qend
        read.mid_seq_end_pos    = read.right_anchor_paf.qstart
        if read.core_seq_start_pos < 0: read.core_seq_start_pos = 0
        if read.core_seq_end_pos > read.full_read_len: read.core_seq_end_pos = read.full_read_len
        read.left_buffer_len = read.left_anchor_paf.qend - read.core_seq_start_pos
        read.right_buffer_len = read.core_seq_end_pos - read.right_anchor_paf.qstart
        if read.left_anchor_paf.strand == '+':
            read.strand = '+'
        else:
            read.strand = '-'

    return

def find_anchor_locations_from_paf(repeat_region: RepeatRegion, anchor_locations_paf_file: string):

    anchor_locations_paf_f = open(anchor_locations_paf_file, 'r')
    read_paf_list = list()
    while 1:
        line = anchor_locations_paf_f.readline()
        if not line: break
        line = line.strip()
        if not line: continue

        col_list = line.strip().split('\t')
        paf = PAF(col_list)

        if len(read_paf_list) == 0 or paf.qname == read_paf_list[0].qname:
            read_paf_list.append(paf)
        else:
            find_anchor_locations_for1read(read_paf_list, repeat_region)
            read_paf_list.clear()
            read_paf_list.append(paf)

    anchor_locations_paf_f.close()
    
    find_anchor_locations_for1read(read_paf_list, repeat_region)

    return

def find_anchor_locations_in_reads(minimap2:string, repeat_region:RepeatRegion, num_cpu:int):

    left_template_seq    = repeat_region.left_anchor_seq
    left_template_name   = 'left_anchor'

    right_template_seq   = repeat_region.right_anchor_seq
    right_template_name  = 'right_anchor'

    template_fasta_file  = os.path.join(repeat_region.temp_out_dir, 'anchors.fasta')
    repeat_region.temp_file_list.append(template_fasta_file)

    template_fasta_fp    = open(template_fasta_file, 'w')
    template_fasta_fp.write('>%s\n' % left_template_name)
    template_fasta_fp.write('%s\n' % left_template_seq)
    template_fasta_fp.write('>%s\n' % right_template_name)
    template_fasta_fp.write('%s\n' % right_template_seq)
    template_fasta_fp.close()

    anchor_locations_paf_file = os.path.join(repeat_region.temp_out_dir, 'anchor_locations.paf')
    repeat_region.temp_file_list.append(anchor_locations_paf_file)

    preset = tk.get_preset_for_minimap2('ont')
    cmd = f'{minimap2} -c -t {num_cpu} -x {preset} {template_fasta_file} {repeat_region.region_fq_file} > {anchor_locations_paf_file} 2> /dev/null'
    tk.eprint(f'NOTICE: running command: {cmd}')
    tk.run_system_cmd(cmd)
    find_anchor_locations_from_paf(repeat_region, anchor_locations_paf_file)
    
    return

def make_core_seq_fastq(repeat_region: RepeatRegion):
    repeat_region.core_seq_fq_file = os.path.join(repeat_region.temp_out_dir, 'core_sequences.fastq')
    repeat_region.temp_file_list.append(repeat_region.core_seq_fq_file)

    repeat_region.mid_seq_fq_file = os.path.join(repeat_region.temp_out_dir, 'middle_sequences.fastq')
    repeat_region.temp_file_list.append(repeat_region.mid_seq_fq_file)

    region_fq_f   = open(repeat_region.region_fq_file, 'r')
    core_seq_fq_f = open(repeat_region.core_seq_fq_file, 'w')
    mid_seq_fq_f  = open(repeat_region.mid_seq_fq_file, 'w')
    while 1:
        line1 = region_fq_f.readline()
        line2 = region_fq_f.readline()
        line3 = region_fq_f.readline()
        line4 = region_fq_f.readline()
        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        read_name = line1.strip()[1:]
        if read_name not in repeat_region.read_dict: continue
        read_sequence = line2.strip()            
        if repeat_region.read_dict[read_name].strand == '-':
            read_sequence = tk.rev_comp(read_sequence)
        core_sequence = read_sequence[repeat_region.read_dict[read_name].core_seq_start_pos:repeat_region.read_dict[read_name].core_seq_end_pos]
        mid_sequence = read_sequence[repeat_region.read_dict[read_name].mid_seq_start_pos:repeat_region.read_dict[read_name].mid_seq_end_pos]

        core_seq_fq_f.write(line1)
        core_seq_fq_f.write(core_sequence + '\n')
        core_seq_fq_f.write(line3)
        core_seq_fq_f.write('0' * len(core_sequence) + '\n')


        mid_seq_fq_f.write(line1)
        mid_seq_fq_f.write(mid_sequence + '\n')
        mid_seq_fq_f.write(line3)
        mid_seq_fq_f.write('0' * len(mid_sequence) + '\n')

    region_fq_f.close()
    core_seq_fq_f.close()
    mid_seq_fq_f.close()
    return

def output_results(repeat_region:RepeatRegion):

    read_count_file =  f'{repeat_region.out_prefix}.repeat_size.txt'
    read_count_f = open(read_count_file, 'w')

    repeat_size_list = []
    for read_name in repeat_region.read_dict:
        repeat_size = repeat_region.read_dict[read_name].round3_repeat_size
        if repeat_size != None:
            repeat_size_list.append(repeat_size)

    if len(repeat_size_list) > 0:
        median_repeat_size = np.median(repeat_size_list)
    else:
        median_repeat_size = 'N.A.'
    
    num_reads = len(repeat_size_list)
    read_count_f.write(f'#predicted_repeat_size={median_repeat_size}\tnum_reads={num_reads}\n')
    for read_name in repeat_region.read_dict:
        repeat_size = repeat_region.read_dict[read_name].round3_repeat_size
        if repeat_size != None:
            read_count_f.write(f'{read_name}\t{repeat_size:.1f}\n')
    read_count_f.close()

def round1_and_round2_estimation(minimap2: string, repeat_region: RepeatRegion, num_cpu: int):

    if len(repeat_region.read_dict) == 0: return

    round1_repeat_size_list = []
    for read_name in repeat_region.read_dict:
        read = repeat_region.read_dict[read_name]
        read.round1_repeat_size = float(read.dist_between_anchors)/len(repeat_region.repeat_unit_seq)
        round1_repeat_size_list.append(read.round1_repeat_size)
    
    template_repeat_size = int(max(round1_repeat_size_list) * 1.5) + 1

    if template_repeat_size < max(round1_repeat_size_list) + 10:
        template_repeat_size = int(max(round1_repeat_size_list) + 10)
    
    round1_fasta_file  = os.path.join(repeat_region.temp_out_dir, 'round1_ref.fasta')
    repeat_region.temp_file_list.append(round1_fasta_file)
    
    round1_fasta_f = open(round1_fasta_file, 'w')
    round1_fasta_f.write(f'>{template_repeat_size}\n')
    round1_fasta_f.write(f'{repeat_region.left_anchor_seq}{repeat_region.repeat_unit_seq * template_repeat_size}\n')
    round1_fasta_f.close()

    round1_paf_file = os.path.join(repeat_region.temp_out_dir, 'round1.paf')
    repeat_region.temp_file_list.append(round1_paf_file)
    minimap2_parameters = ' -x map-ont -B 9 '
    cmd = f'{minimap2} -c -t {num_cpu} {minimap2_parameters} {round1_fasta_file} {repeat_region.core_seq_fq_file} > {round1_paf_file} 2> /dev/null'
    
    tk.eprint(f'NOTICE: running command: {cmd}')
    tk.run_system_cmd(cmd)
    round1_paf_f = open(round1_paf_file, 'r')
    lines = list(round1_paf_f)
    round1_paf_f.close()

    for line in lines:
        col_list = line.strip().split('\t')
        paf = PAF(col_list)
        if paf.tstart <= len(repeat_region.left_anchor_seq) and paf.tend >= len(repeat_region.left_anchor_seq):
            read_name = paf.qname
            round2_repeat_size = float(paf.tend - len(repeat_region.left_anchor_seq))/len(repeat_region.repeat_unit_seq)
            repeat_region.read_dict[read_name].round2_repeat_size = round2_repeat_size
    
    return

def quantify1repeat_from_bam(input_args, in_bam_file, ref_fasta_dict, repeat_region:RepeatRegion):

    temp_out_dir = f'{input_args.out_prefix}.NanoRepeat_temp_dir.{repeat_region.to_unique_id()}'
    os.makedirs(temp_out_dir, exist_ok=True)

    repeat_region.temp_out_dir = temp_out_dir
    repeat_region.out_prefix = f'{input_args.out_prefix}.{repeat_region.to_unique_id()}'
    repeat_region.anchor_len = input_args.anchor_len

    # extract reads from bam file
    repeat_region.region_fq_file = os.path.join(temp_out_dir, f'{repeat_region.to_unique_id()}.fastq')
    cmd = f'{input_args.samtools} view -h {in_bam_file} {repeat_region.to_invertal(flank_dist=10)} | {input_args.samtools} fastq - > {repeat_region.region_fq_file}'
    tk.run_system_cmd(cmd)

    # extract ref sequence
    extract_ref_sequence(ref_fasta_dict, repeat_region)
    
    tk.eprint('NOTICE: step 1: finding anchor location in reads')
    find_anchor_locations_in_reads(input_args.minimap2, repeat_region, input_args.num_cpu)
    tk.eprint('NOTICE: step 1 finished')

    # make core sequence fastq
    make_core_seq_fastq(repeat_region)

    tk.eprint('NOTICE: step 2: round 1 and round 2 estimation')
    round1_and_round2_estimation(input_args.minimap2, repeat_region, input_args.num_cpu)
    tk.eprint('NOTICE: step 2 finished')

    tk.eprint('NOTICE: step 3: round 3 estimation')
    round3_estimation(input_args.minimap2, repeat_region, input_args.num_cpu)
    tk.eprint('NOTICE: step 3 finished')

    split_allele_using_gmm(repeat_region, input_args.ploidy)
    
    tk.eprint('NOTICE: program finished.')
    shutil.rmtree(repeat_region.temp_out_dir)

    return

def split_allele_using_gmm(repeat_region, ploidy):

    if ploidy < 1:
        tk.eprint('ERROR! ploidy must be >= 1 !')
        sys.exit()

    read_repeat_count_dict = dict()
    for read_name in repeat_region.read_dict:
        if repeat_region.read_dict[read_name].round3_repeat_size != None:
            read_repeat_count_dict[read_name] = repeat_region.read_dict[read_name].round3_repeat_size

    in_fastq_file = repeat_region.region_fq_file
    out_prefix = repeat_region.out_prefix

    in_fastq_filename = os.path.split(in_fastq_file)[1]
    proba_cutoff = 0.95
    cov_type = 'tied'
    
    min_repeat_count_cutoff, max_repeat_count_cutoff = analysis_outlier(read_repeat_count_dict)
    readname_list = list()
    read_repeat_count_list = list()
    for readname in read_repeat_count_dict:
        repeat_count = read_repeat_count_dict[readname]
        if repeat_count < min_repeat_count_cutoff or repeat_count > max_repeat_count_cutoff: continue 
        readname_list.append(readname)
        read_repeat_count_list.append(repeat_count)

    read_repeat_count_array = np.array(read_repeat_count_list)
    num_data_points = len(read_repeat_count_list)
    read_repeat_count_array = read_repeat_count_array.reshape(num_data_points, 1)

    if num_data_points < 1:
        sys.stderr.write('ERROR! no enough reads!')
        sys.stderr.write('ERROR! number of reads passed filtering: %d' % num_data_points)

    best_n_components = chose_best_num_components (read_repeat_count_array, ploidy, proba_cutoff, cov_type)
    final_gmm = GaussianMixture(n_components=best_n_components, covariance_type=cov_type, n_init = 10).fit(read_repeat_count_array)
    old_read_label_list = list(final_gmm.predict(read_repeat_count_array))
    proba2darray = final_gmm.predict_proba(read_repeat_count_array)
    old_cluster_mean_list = list(final_gmm.means_)    

    read_label_list, old_label_to_new_label_dict, new_label_to_old_label_dict = sort_label_by_cluster_mean(old_read_label_list, old_cluster_mean_list)

    read_label_dict = dict()
    read_proba_dict = dict()
    qc_failed_readname_set = set()
    all_info_dict = dict()
    
    each_allele_repeat_count_2d_list = [0] * best_n_components
    for i in range(0, len(each_allele_repeat_count_2d_list)):
        each_allele_repeat_count_2d_list[i] = list()

    for i in range(0, len(readname_list)):
        readname = readname_list[i]
        repeat_count = read_repeat_count_list[i]
        read_label = read_label_list[i]
        read_label_dict[readname] = read_label
        proba_array = proba2darray[i]
        max_prob = max(proba_array)
        read_proba_dict[readname] = max_prob
        info = (repeat_count, read_label, max_prob)
        all_info_dict[readname] = info
        each_allele_repeat_count_2d_list[read_label].append(repeat_count)
        
    for read_label in range(0, len(each_allele_repeat_count_2d_list)):
        allele_repeat_count_list = each_allele_repeat_count_2d_list[read_label]
        old_label = new_label_to_old_label_dict[read_label]
        gmm_average_repeat_number = old_cluster_mean_list[old_label]
        std = np.std(allele_repeat_count_list)
        max_repeat_number = gmm_average_repeat_number + 3 * std
        min_repeat_number = gmm_average_repeat_number - 3 * std

        for readname in all_info_dict:
            repeat_count, label, max_prob = all_info_dict[readname]
            if max_prob < proba_cutoff: qc_failed_readname_set.add(readname)
            if label == read_label:
                if repeat_count > max_repeat_number or repeat_count < min_repeat_number:
                    qc_failed_readname_set.add(readname)

    qc_passed_each_allele_repeat_count_2d_list = [0] * best_n_components
    for i in range(0, len(qc_passed_each_allele_repeat_count_2d_list)):
        qc_passed_each_allele_repeat_count_2d_list[i] = list()
    for readname in all_info_dict:
        repeat_count, label, max_prob = all_info_dict[readname]
        if readname not in qc_failed_readname_set:
            qc_passed_each_allele_repeat_count_2d_list[label].append(repeat_count)

    out_summary_file = out_prefix + '.summary.txt'
    out_detail_file  = out_prefix + '.details.txt'
    hist_figure_file = out_prefix + '.hist.png'

    out_allele_fastq_file_list = list()
    for label in range(0, best_n_components):
        allele_id = label + 1
        out_allele_fastq_file = out_prefix + '.allele%d.fastq' % (allele_id)
        out_allele_fastq_file_list.append(out_allele_fastq_file)

    out_allele_fastq_fp_list = list()
    for i in range(0, len(out_allele_fastq_file_list)):
        out_allele_fastq_fp = open(out_allele_fastq_file_list[i], 'w')
        out_allele_fastq_fp_list.append(out_allele_fastq_fp)

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
        if readname not in read_label_dict: continue
        if readname  in qc_failed_readname_set: continue
    
        label = read_label_dict[readname]
        out_allele_fastq_fp_list[label].write(line1 + line2 + line3 + line4)

    in_fastq_fp.close()

    for i in range(0, len(out_allele_fastq_fp_list)):
        out_allele_fastq_fp_list[i].close()

    predicted_repeat_count_list = list()
    meta_info_list = list()
    out_summray_fp = open(out_summary_file, 'w')
    out_detail_fp  = open(out_detail_file, 'w')
    out_summray_fp.write(f'Repeat_Region={repeat_region.to_unique_id()}\n')
    for read_label in range(0, len(qc_passed_each_allele_repeat_count_2d_list)):
        allele_repeat_count_list = qc_passed_each_allele_repeat_count_2d_list[read_label]
        allele_id = read_label + 1
        meta_info = Metainfo().init_from_repeat_size_list(allele_repeat_count_list, allele_id)
        meta_info_list.append(meta_info)
        predicted_repeat_count_list.append(meta_info.predicted_repeat_size)
        out_summray_fp.write(f'{meta_info.output()}\n')
    
  
    out_info_list = list()
    for readname in all_info_dict:
        if readname in qc_failed_readname_set: continue
        repeat_count, label, max_prob = all_info_dict[readname]
        allele_id = label + 1
        out_info = (readname, repeat_count, allele_id, max_prob)
        out_info_list.append(out_info)

    out_info_list.sort(key = lambda x:x[1])

    out_detail_fp.write('#readname\trepeat_size\tallele_id\tprobablity\tconfidence_of_allele_assignment\n')
    for i in range(0, len(out_info_list)):
        readname, repeat_count, allele_id, max_prob = out_info_list[i]
        if max_prob < proba_cutoff: 
            confidence = 'LOW'
        else:
            confidence = 'HIGH'
        out_detail_fp.write(f'{readname}\t{repeat_count:.1f}\t{allele_id}\t{max_prob:.2f}\t{confidence}\n')
        #out_detail_fp.write('%s\t%d\t%d\t%f\n' % (readname, repeat_count, allele_id, max_prob))

    out_summray_fp.close()
    out_detail_fp.close()
    plot_repeat_counts(qc_passed_each_allele_repeat_count_2d_list, predicted_repeat_count_list, hist_figure_file)

    return


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

def chose_best_num_components (read_repeat_count_array, ploidy, proba_cutoff, cov_type):

    num_useful_data_points_list = list()
    num_useful_data_points_list.append(0)
    
    max_num_components = ploidy

    bic_list = list()
    bic_list.append(0)

 
    for n in range(1, max_num_components+1):
        gmm = GaussianMixture(n_components=n, covariance_type=cov_type, n_init=10).fit(read_repeat_count_array)
        bic = gmm.bic(read_repeat_count_array)
        bic_list.append(bic)

    min_bic = 1e99
    min_bic_n_components = 0
    for i in range(1, len(bic_list)):
        if bic_list[i] < min_bic:
            min_bic = bic_list[i]
            min_bic_n_components = i

    return min_bic_n_components

def analysis_outlier(read_repeat_count_dict):

    read_repeat_count_list = list()
    for readname in read_repeat_count_dict:
        repeat_count = read_repeat_count_dict[readname]
        read_repeat_count_list.append(repeat_count)

    try:
        mean = np.mean(read_repeat_count_list)
    except:
        print(read_repeat_count_list)
        sys.exit(1)
    std = np.std(read_repeat_count_list)

    min_repeat_count_cutoff = mean - 3 * std
    if min_repeat_count_cutoff < 0: min_repeat_count_cutoff = 0
    max_repeat_count_cutoff = mean + 3 * std
 
    return min_repeat_count_cutoff, max_repeat_count_cutoff

def plot_repeat_counts(each_allele_repeat_count_2d_list, predicted_repeat_count_list, out_file):

    plt.figure(figsize=(6, 4))

    all_data_list = list()
    for x in each_allele_repeat_count_2d_list:
        all_data_list += x

    if len(all_data_list) == 0:
        return

    xmin = int(min(all_data_list))
    xmax = int(max(all_data_list))+1

    margin = int((xmax-xmin+1)/2)
    xmin -= margin
    xmax += margin
    if xmin < 0: xmin = 0

    if xmax - xmin < 200:
        b = range(xmin - 1, xmax + 2)
    else:
        b = range(xmin - 1, xmax + 2, int((xmax - xmin)/200))

    for i in range(0, len(each_allele_repeat_count_2d_list)):
        x = each_allele_repeat_count_2d_list[i]
        plt.hist(x, bins = b)

    for repeat_count in predicted_repeat_count_list:
        plt.axvline(x=repeat_count, color = 'grey', linestyle = ':')

    plt.title('Repeat size distribution')
    plt.xlabel('repeat size')
    plt.ylabel('number of reads')
    plt.savefig(out_file, dpi=300)
    plt.close('all')
    
    return

def get_max_read_length_from_fastq(fq_file):

    max_read_length = 0
    fq_f = open(fq_file, 'r')
    while 1:
        line = fq_f.readline()
        if not line: break
        if len(line)-1 > max_read_length:
            max_read_length = len(line)-1
    fq_f.close()

    return max_read_length

def round3_estimation_for1read(repeat_region: RepeatRegion, read_paf_list):

    if len(read_paf_list) == 0: return

    read_paf_list.sort(key = lambda paf:paf.align_score, reverse = True)
    best_paf_repeat_size_list = []
    read_name = read_paf_list[0].qname
    for i in range(0, len(read_paf_list)):
        if read_paf_list[i].align_score < read_paf_list[0].align_score: break
        if read_paf_list[i].tstart < repeat_region.left_anchor_len and read_paf_list[i].tlen - read_paf_list[i].tend < repeat_region.right_anchor_len and read_paf_list[i].align_score == read_paf_list[0].align_score:
            best_paf_repeat_size_list.append(int(read_paf_list[i].tname))

    if len(best_paf_repeat_size_list) > 0:
        repeat_region.read_dict[read_name].round3_repeat_size = np.mean(best_paf_repeat_size_list)
    else:
        repeat_region.read_dict[read_name].round3_repeat_size = repeat_region.read_dict[read_name].round2_repeat_size
    return
                
def round3_estimation_from_paf(repeat_region: RepeatRegion, round3_paf_file):

    round3_paf_f = open(round3_paf_file, 'r')
    lines = list(round3_paf_f)
    round3_paf_f.close()

    read_paf_list = []
    for line in lines:
        col_list = line.strip().split('\t')
        paf = PAF(col_list)
        if paf.qname not in repeat_region.read_dict: continue
        if repeat_region.read_dict[paf.qname].round2_repeat_size == None: continue
        if len(read_paf_list) == 0 or paf.qname == read_paf_list[0].qname:
            read_paf_list.append(paf)
        else:
            round3_estimation_for1read(repeat_region, read_paf_list)
            read_paf_list.clear()
            read_paf_list.append(paf)

    round3_estimation_for1read(repeat_region, read_paf_list)
    
    return 


def round3_estimation(minimap2:string, repeat_region:RepeatRegion, num_cpu:int):
    
    estimated_repeat_size_list = []
    for read_name in repeat_region.read_dict:
        read = repeat_region.read_dict[read_name]
        if read.round2_repeat_size != None:
            estimated_repeat_size_list.append(read.round2_repeat_size)
    if len(estimated_repeat_size_list) == 0: return
    max_template_repeat_size = int(max(estimated_repeat_size_list) * 1.5) + 1
    min_template_repeat_size = int(min(estimated_repeat_size_list) / 2)

    if min_template_repeat_size > min(estimated_repeat_size_list) - 10:
        min_template_repeat_size = int(min(estimated_repeat_size_list) - 10)
    if max_template_repeat_size < max(estimated_repeat_size_list) + 10:
        max_template_repeat_size = int(max(estimated_repeat_size_list) + 10) + 1
    
    if min_template_repeat_size < 0: min_template_repeat_size = 0
    
    template_fasta_file  = os.path.join(repeat_region.temp_out_dir, 'round3_ref.fasta')
    repeat_region.temp_file_list.append(template_fasta_file)
    template_fasta_fp    = open(template_fasta_file, 'w')
    for repeat_size in range(min_template_repeat_size, max_template_repeat_size+1):
        template_seq = repeat_region.left_anchor_seq + repeat_region.repeat_unit_seq * repeat_size + repeat_region.right_anchor_seq
        template_fasta_fp.write('>%s\n' % repeat_size)
        template_fasta_fp.write('%s\n' % template_seq)

    template_fasta_fp.close()

    round3_paf_file = os.path.join(repeat_region.temp_out_dir, 'round3.paf')
    repeat_region.temp_file_list.append(round3_paf_file)

    minimap2_parameters = ' -x map-ont -B 9 '
    cmd = f'{minimap2} {minimap2_parameters} -N 100 -c --eqx -t {num_cpu} {template_fasta_file} {repeat_region.core_seq_fq_file} > {round3_paf_file} 2> /dev/null'
    tk.run_system_cmd(cmd)

    round3_estimation_from_paf(repeat_region, round3_paf_file)
    
    return

class Metainfo:
    def __init__(self, allele_id = -1, num_reads = 0, predicted_repeat_size = -1, min_repeat_size = -1, max_repeat_size = -1):
        self.allele_id = allele_id
        self.num_reads = num_reads
        self.predicted_repeat_size = predicted_repeat_size
        self.min_repeat_size = min_repeat_size
        self.max_repeat_size = max_repeat_size
    
    def init_from_repeat_size_list(self, repeat_size_list, allele_id):
        
        self.allele_id = allele_id
        self.num_reads = len(repeat_size_list)
        if len(repeat_size_list) == 0:
            self.predicted_repeat_size = -1
            self.min_repeat_size = -1 
            self.max_repeat_size = -1
            return self
        
        self.min_repeat_size = min(repeat_size_list)
        self.max_repeat_size = max(repeat_size_list)
        self.predicted_repeat_size = int(np.median(repeat_size_list) + 0.5)
        return self

    def output(self):
        return 'allele_id=%d\tnum_reads=%d\tpredicted_repeat_size=%d\tmin_repeat_size=%d\tmax_repeat_size=%d' % (self.allele_id, self.num_reads, self.predicted_repeat_size, self.min_repeat_size, self.max_repeat_size)

def nanoRepeat_bam (input_args, in_bam_file):
    
    tk.eprint(f'NOTICE: reading repeat region file: {input_args.repeat_region_bed}')
    repeat_region_list = read_repeat_region_file(input_args.repeat_region_bed)

    tk.eprint(f'NOTICE: reading referece fasta file: {input_args.ref_fasta}')
    ref_fasta_dict = tk.fasta_file2dict(input_args.ref_fasta)
    
    for repeat_region in repeat_region_list:
        tk.eprint(f'NOTICE: quantifying repeat: {repeat_region.to_unique_id()}')
        quantify1repeat_from_bam(input_args, in_bam_file, ref_fasta_dict, repeat_region)

    return

