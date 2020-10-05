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
import numpy as np
import math
import argparse
import subprocess

from numpy.f2py.auxfuncs import containsmodule
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip

import tk


tab  = '\t' 
endl = '\n'

def parse_user_arguments():

    parser = argparse.ArgumentParser(description='A tools for short tandem repeat detection from long-read amplicon sequencing data')
    ### required arguments ###
    parser.add_argument('--in_fq',              required = True, metavar = 'PATH',   type = str, help = 'path to input fastq file')
    parser.add_argument('--platform',           required = True, metavar = 'STRING', type = str, help = 'three valid values: ont, pacbio, consensus')
    parser.add_argument('--ref_fasta',          required = True, metavar = 'PATH',   type = str, help = 'path to reference genome sequence in FASTA format')
    parser.add_argument('--repeat_region',      required = True, metavar = 'chr:start-end', type = str, help = 'repeat region in the reference genome (e.g. chr4:3074876-3074939, coordinates are 1-based)')
    parser.add_argument('--repeat_unit',        required = True, metavar = 'STRING', type = str, help = 'sequence of the repeat unit (e.g. CAG)')
    parser.add_argument('--max_repeat_size',    required = True, metavar = 'INT',    type = int, help = 'maximum possible number of the repeat unit')
    parser.add_argument('--out_dir',            required = True, metavar = 'PATH',   type = str, help = 'path to the output directory')
   
    ### optional arguments ###
    parser.add_argument('--num_threads',        required = False, metavar = 'INT',   type = int, default = 1,  help ='number of threads used by minimap2 (default: 1)')
    parser.add_argument('--minimap2',           required = False, metavar = 'PATH',  type = str, default = '', help ='path to minimap2 (default: using environment default)')
    parser.add_argument('--fixed_cutoff_value', required = False, metavar = 'INT',   type = int, default = -1, help ='split alleles using this fixed_cutoff_value (if set, ploidy will be set to 2). reads with repeat size >= fixed_cutoff_value will be assigned to the second allele')
    parser.add_argument('--ploidy',             required = False, metavar = 'INT',   type = int, default = 2,  help ='ploidy of the sample (default: 2)')
    parser.add_argument('--anchor_len',         required = False, metavar = 'INT',   type = int, default = 1000, help ='length of up/downstream sequence to help identify the repeat region (default: 1000 bp, increase this value if the 1000 bp up/downstream sequences are also repeat)')
    parser.add_argument('--version', action='version', version='%(prog)s 0.3.0')

    input_args = parser.parse_args()
    if input_args.ploidy <= 0: 
        tk.eprint('ERROR! ploidy must be > 0')
        sys.exit()

    input_args.platform = input_args.platform.strip('`')
    input_args.platform = input_args.platform.strip('\'')
    input_args.platform = input_args.platform.strip('\"')

    valid_platform_list = ['ont', 'pacbio', 'consensus']
    input_args.platform = input_args.platform.lower()

    if input_args.platform not in valid_platform_list:
        tk.eprint('ERROR! --platform should be one of the three valid values: ont, pacbio, consensus !\n')
        sys.exit()

    tk.check_input_file_exists(input_args.in_fq)
    tk.check_input_file_exists(input_args.ref_fasta)

    input_args.in_fq     = os.path.abspath(input_args.in_fq)
    input_args.ref_fasta = os.path.abspath(input_args.ref_fasta)
    input_args.out_dir   = os.path.abspath(input_args.out_dir)
    if input_args.minimap2 == '': 
        input_args.minimap2 = tk.find_executable_path('minimap2')
        if not input_args.minimap2:
            tk.eprint('ERROR! minimap2 was not found! Please supply the path to --minimap2')
            sys.exit()
        else:
            if os.path.exists(input_args.minimap2):
                tk.eprint('NOTICE: found path to minimap2: %s' % input_args.minimap2)
            else:
                tk.eprint('ERROR! minimap2 was not found! Please supply the path to --minimap2')

    tk.check_input_file_exists(input_args.minimap2)

    return input_args

def main():

    input_args = parse_user_arguments()   
    ampRepeat (input_args)

    return

def ampRepeat (input_args):

    minimap2            = input_args.minimap2
    ref_fasta_file      = input_args.ref_fasta
    repeat_unit         = input_args.repeat_unit
    repeat_region       = input_args.repeat_region
    in_fastq_file       = input_args.in_fq
    platform            = input_args.platform
    max_repeat_size     = input_args.max_repeat_size
    num_threads         = input_args.num_threads
    out_dir             = input_args.out_dir
    max_anchor_len      = input_args.anchor_len
    fixed_cutoff_value  = input_args.fixed_cutoff_value
    ploidy              = input_args.ploidy


    tk.eprint('NOTICE: input fastq file is: %s' % in_fastq_file)
    tk.eprint('NOTICE: ref fasta file is: %s' % ref_fasta_file)
    tk.eprint('NOTICE: output dir is: %s' % out_dir)
    tk.create_dir(out_dir)
    in_fastq_prefix = os.path.splitext(os.path.split(in_fastq_file)[1])[0]
    temp_out_dir = os.path.join(out_dir, '%s.AmpRepeat_temp' % in_fastq_prefix)
    tk.create_dir(temp_out_dir)

    repeat_chr, repeat_start, repeat_end = analysis_repeat_region(repeat_region)

    tk.eprint('NOTICE: intitial estimation of repeat size...')
    first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set = initial_estimation(in_fastq_file, ref_fasta_file, repeat_unit, repeat_chr, repeat_start, repeat_end, max_repeat_size, minimap2, platform, num_threads, max_anchor_len, temp_out_dir)
    
    fine_tuned_repeat_count_dict = fine_tune_read_count(first_round_repeat_count_dict, potential_repeat_region_dict, in_fastq_file, ref_fasta_file, repeat_unit, repeat_chr, repeat_start, repeat_end, max_repeat_size, minimap2, platform, num_threads, max_anchor_len, temp_out_dir)

    for readname in bad_reads_set:
        fine_tuned_repeat_count_dict.pop(readname, None)

    tk.eprint('NOTICE: assigned repeat size to %d reads.' % len(fine_tuned_repeat_count_dict))
    if fixed_cutoff_value >= 0:
        split_allele_using_fixed_cutoff_value (fixed_cutoff_value, fine_tuned_repeat_count_dict, in_fastq_file, out_dir)
    else:
        split_allele_using_gmm(ploidy, fine_tuned_repeat_count_dict, in_fastq_file, out_dir)
    
    cmd = f'rm -r {temp_out_dir}'
    tk.run_system_cmd(cmd)
    tk.eprint('NOTICE: program finished. output files are here: %s' % out_dir)

    return

def split_allele_using_fixed_cutoff_value (fixed_cutoff_value, read_repeat_count_dict, in_fastq_file, out_dir):

    in_fastq_prefix = os.path.splitext(os.path.split(in_fastq_file)[1])[0]
    in_fastq_filename = os.path.split(in_fastq_file)[1]
    out_prefix = '%s.fixed_cutoff_%d' % (in_fastq_prefix, fixed_cutoff_value)

    out_allele1_fastq_file = os.path.join(out_dir, '%s.allele1.fastq' % (out_prefix))
    out_allele2_fastq_file = os.path.join(out_dir, '%s.allele2.fastq' % (out_prefix))
    out_summray_file       = os.path.join(out_dir, '%s.summary.txt' % (out_prefix))
    hist_figure_file       = os.path.join(out_dir, '%s.hist.png' % (out_prefix))

    allele1_repeat_count_list = list()
    allele2_repeat_count_list = list()

    in_fastq_fp = tk.gzopen(in_fastq_file)
    out_allele1_fastq_fp = open(out_allele1_fastq_file, 'w')
    out_allele2_fastq_fp = open(out_allele2_fastq_file, 'w')
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
        if readname not in read_repeat_count_dict: continue

        repeat_count = read_repeat_count_dict[readname]
        if repeat_count < fixed_cutoff_value:
            allele1_repeat_count_list.append(repeat_count)
            out_allele1_fastq_fp.write(line1 + line2 + line3 + line4)
        else:
            allele2_repeat_count_list.append(repeat_count)
            out_allele2_fastq_fp.write(line1 + line2 + line3 + line4)

    in_fastq_fp.close()
    out_allele1_fastq_fp.close()
    out_allele2_fastq_fp.close()

    read_repeat_count_list = list()
    for readname, repeat_size in read_repeat_count_dict.items():
        if repeat_size < fixed_cutoff_value:
            allele_id = 1
        else:
            allele_id = 2
        tup = (readname, repeat_size, allele_id)
        read_repeat_count_list.append(tup)
    read_repeat_count_list.sort(key = lambda x:x[1])

    metainfo1 = Metainfo().init_from_repeat_size_list(allele1_repeat_count_list, 1)
    metainfo2 = Metainfo().init_from_repeat_size_list(allele2_repeat_count_list, 2)

    out_summray_fp = open(out_summray_file, 'w')
    header = f'##{in_fastq_filename}\tfixed_cutoff={fixed_cutoff_value}\t{metainfo1.output()}\t{metainfo2.output()}\n'
    out_summray_fp.write(header)
    out_summray_fp.write('#readname\trepeat_size\tallele\n')
    for tup in read_repeat_count_list:
        out_summray_fp.write(f'{tup[0]}\t{int(tup[1])}\t{tup[2]}\n')
    out_summray_fp.close()

    each_allele_repeat_count_2d_list = [allele1_repeat_count_list, allele2_repeat_count_list]
    predicted_repeat_count_list      = [metainfo1.predicted_repeat_size, metainfo2.predicted_repeat_size]

    plot_repeat_counts(each_allele_repeat_count_2d_list, predicted_repeat_count_list, hist_figure_file)

    return

def split_allele_using_gmm(ploidy, read_repeat_count_dict, in_fastq_file, out_dir):

    if ploidy < 1:
        tk.eprint('ERROR! ploidy must be >= 1 !')
        sys.exit()

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
        sys.exit()

    best_n_components = chose_best_num_components (read_repeat_count_array, ploidy, proba_cutoff, cov_type)
    #best_n_components = ploidy
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


    in_fastq_prefix = os.path.splitext(os.path.split(in_fastq_file)[1])[0]
 
    out_prefix = os.path.join(out_dir, '%s.GMM' % (in_fastq_prefix))

    out_summray_file = out_prefix + '.summary.txt'
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
    out_summray_fp = open(out_summray_file, 'w')
    header = f'##{in_fastq_filename}\tGMM'
    for read_label in range(0, len(qc_passed_each_allele_repeat_count_2d_list)):
        allele_repeat_count_list = qc_passed_each_allele_repeat_count_2d_list[read_label]
        allele_id = read_label + 1
        meta_info = Metainfo().init_from_repeat_size_list(allele_repeat_count_list, allele_id)
        meta_info_list.append(meta_info)
        predicted_repeat_count_list.append(meta_info.predicted_repeat_size)
        header += f'\t{meta_info.output()}'
    
    out_summray_fp.write(header + '\n')
  
    out_info_list = list()
    for readname in all_info_dict:
        if readname in qc_failed_readname_set: continue
        repeat_count, label, max_prob = all_info_dict[readname]
        allele_id = label + 1
        out_info = (readname, repeat_count, allele_id, max_prob)
        out_info_list.append(out_info)

    out_info_list.sort(key = lambda x:x[1])

    out_summray_fp.write('#readname\trepeat_size\tallele\tprobablity\n')
    for i in range(0, len(out_info_list)):
        readname, repeat_count, allele_id, max_prob = out_info_list[i]
        if max_prob < proba_cutoff: continue
        out_summray_fp.write('%s\t%d\t%d\t%f\n' % (readname, repeat_count, allele_id, max_prob))

    out_summray_fp.close()

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

    mean = np.mean(read_repeat_count_list)
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

def fine_tune_read_count(first_round_repeat_count_dict, potential_repeat_region_dict, in_fastq_file, ref_fasta_file, repeat_unit, repeat_chr, repeat_start, repeat_end, max_repeat_size, minimap2, platform, num_threads, max_anchor_len, out_dir):

    first_round_max_repeat_size = 0
    for readname, value in first_round_repeat_count_dict.items():
        min_repeat_size, max_repeat_size = value
        if max_repeat_size > first_round_max_repeat_size:
            first_round_max_repeat_size = max_repeat_size
    
    tk.eprint(f'NOTICE: max repeat size from first round estimation is {first_round_max_repeat_size}')

    tk.eprint('NOTICE: fine-tuning repeat size')
    fine_tuned_repeat_count_dict = dict()

    repeat_chrom_seq  = tk.read_one_chr_from_fasta_file(ref_fasta_file, repeat_chr)
    left_anchor_seq, right_anchor_seq = extract_anchor_sequence(repeat_chrom_seq, repeat_start, repeat_end, max_anchor_len)
    left_anchor_len = len(left_anchor_seq)
    right_anchor_len = len(right_anchor_seq)
    repeat_unit_size = len(repeat_unit)
    fastq_chunk_list = split_input_fastq_file_by_repeat_size(first_round_repeat_count_dict, potential_repeat_region_dict, in_fastq_file, first_round_max_repeat_size, out_dir)
    
    preset = tk.get_preset_for_minimap2(platform)

    tk.eprint('NOTICE: second round alignment')
    chunk_ref_fasta_file_list = list()

    for repeat_size in range(0, first_round_max_repeat_size+10):
        chunk_ref_fasta_file = os.path.join(out_dir, f'{repeat_size}.ref.fasta')
        chunk_ref_fasta_fp = open(chunk_ref_fasta_file, 'w')
        template_seq = left_anchor_seq + repeat_unit * repeat_size + right_anchor_seq
        chunk_ref_fasta_fp.write(f'>{repeat_size}\n')
        chunk_ref_fasta_fp.write(template_seq + '\n')
        chunk_ref_fasta_fp.close()
        chunk_ref_fasta_file_list.append(chunk_ref_fasta_file)
    
    for fastq_chunk in fastq_chunk_list:
        
        max_size = fastq_chunk.max_size
        min_size = fastq_chunk.min_size
    
        tk.eprint('NOTICE: current chunk is: %d-%d' % (min_size, max_size))
        chunk_paf_file = os.path.join(out_dir, f'{min_size}-{max_size}.paf')
    
        chunk_paf_fp = open(chunk_paf_file, 'w')
        chunk_paf_fp.close()
        
        for repeat_size in range(min_size, max_size):
            if repeat_size >= len(chunk_ref_fasta_file_list): continue
            chunk_ref_fasta_file = chunk_ref_fasta_file_list[repeat_size]
            if repeat_size * repeat_unit_size < 100:
                additional_para = tk.minimap2_mid_para
            else:
                additional_para = ''
            cmd = f'{minimap2} {additional_para} -A 2 -B 5 -N 1 -c --eqx --cs -t {num_threads} -x {preset} {chunk_ref_fasta_file} {fastq_chunk.fn} >> {chunk_paf_file} 2> /dev/null'
            tk.run_system_cmd(cmd)
        
        second_round_estimation_from_paf(chunk_paf_file, left_anchor_len, right_anchor_len, repeat_unit, fine_tuned_repeat_count_dict)
        tk.rm(fastq_chunk.fn)
    
    return fine_tuned_repeat_count_dict

def second_round_estimation_from_paf(second_round_paf_file, left_anchor_len, right_anchor_len, repeat_unit, second_round_repeat_count_dict):

    repeat_unit_size = len(repeat_unit)
    second_round_paf_fp = open(second_round_paf_file, 'r')
    read_align_score_dict = dict()
    while 1:
        line = second_round_paf_fp.readline()
        if not line: break
        line = line.strip()
        if not line: continue

        col_list = line.strip().split('\t')
        paf = tk.PAF(col_list)
        repeat_size = int(paf.tname)
        a = left_anchor_len - 7
        if a < 0: a = 0
        b = left_anchor_len + repeat_unit_size * repeat_size + 7
        if b > paf.tlen: b = paf.tlen
        repeat_region_align_score = tk.target_region_alignment_stats_from_cigar(paf.cigar, paf.tstart, paf.tend, a, b).score
    
        tup = (repeat_size, paf.qlen, repeat_region_align_score, paf.align_score)
        if paf.qname not in read_align_score_dict:
            read_align_score_dict[paf.qname] = list()
        read_align_score_dict[paf.qname].append(tup)    

    second_round_paf_fp.close()

    for readname in read_align_score_dict:
        tup_list = read_align_score_dict[readname]
        tup_list.sort(key = lambda x:x[2], reverse = True)
        max_score = tup_list[0][2]
        max_score_repeat_size_list = list()
        for tup in tup_list:
            if tup[2] == max_score:
                max_score_repeat_size_list.append(tup[0])
            else:
                break
        round2_repeat_size = np.mean(max_score_repeat_size_list)
        second_round_repeat_count_dict[readname] = round2_repeat_size
    
    return

def split_input_fastq_file_by_repeat_size(first_round_max_repeat_count_dict, potential_repeat_region_dict, in_fastq_file, max_repeat_size, out_dir):

    tk.eprint('NOTICE: splitting fastq files')
    bin_size = int(max_repeat_size / 200) + 1
    fastq_chunk_list = list()
    repeat_size2chunk_table = [-1] * (max_repeat_size+bin_size+1)
    for min_size in range(0, max_repeat_size+1, bin_size):
        max_size = min_size + bin_size
        fn = os.path.join(out_dir, f'{min_size}-{max_size}.fastq')
        fp = open(fn, 'w')
        fastq_chunk = FastqChunk(min_size, max_size, fn, fp)
        fastq_chunk_list.append(fastq_chunk)

        for i in range(min_size, max_size):
            repeat_size2chunk_table[i] = len(fastq_chunk_list) - 1

    in_fastq_fp = tk.gzopen(in_fastq_file, 'r')
    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()
        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        if line1[0] != '@' or len(line2)!= len(line4):
            tk.eprint(f'ERROR! found corrupted reads in input fastq file: {in_fastq_file}')
            sys.exit()
        
        readname = line1[1:].split()[0]
        if readname not in first_round_max_repeat_count_dict: continue
        min_size, max_size = first_round_max_repeat_count_dict[readname]
        chunk_id1 = repeat_size2chunk_table[min_size]
        chunk_id2 = repeat_size2chunk_table[max_size]

        if chunk_id1 < 0: chunk_id1 = 0
        if chunk_id2 < 0: chunk_id2 = len(fastq_chunk_list)-1
    
        if readname not in potential_repeat_region_dict:
            for x in range(chunk_id1, chunk_id2+1):
                fastq_chunk_list[x].fp.write(line1 + line2 + line3 + line4)
        else:
            start, end = potential_repeat_region_dict[readname]
            seq  = line2[start:end] + '\n'
            qual = line4[start:end] + '\n'
            out = f'@{readname} {start}-{end}\n' + seq + line3 + qual
            for x in range(chunk_id1, chunk_id2+1):
                fastq_chunk_list[x].fp.write(out)
     
    in_fastq_fp.close()
    for fastq_chunk in fastq_chunk_list:
        fastq_chunk.fp.close()

    return fastq_chunk_list

def initial_estimation(in_fastq_file, ref_fasta_file, repeat_unit, repeat_chr, repeat_start, repeat_end, max_repeat_size, minimap2, platform, num_threads, max_anchor_len, out_dir):
    
    repeat_chrom_seq  = tk.read_one_chr_from_fasta_file(ref_fasta_file, repeat_chr)
    left_anchor_seq, right_anchor_seq = extract_anchor_sequence(repeat_chrom_seq, repeat_start, repeat_end, max_anchor_len)
    
    ## make template fasta file ##
    template_repeat_size = max_repeat_size + 10
    template_fasta_file = os.path.join(out_dir, 'templates.fasta')
   
    left_template_seq = left_anchor_seq + repeat_unit * template_repeat_size
    left_template_name = 'left_anchor_%d_%d' % (len(left_anchor_seq), template_repeat_size)
    
    right_template_seq = repeat_unit * template_repeat_size + right_anchor_seq
    right_template_seq = tk.rev_comp(right_template_seq)
    right_template_name = 'right_anchor_%d_%d_revc' % (len(right_anchor_seq), template_repeat_size)
  
    template_fasta_fp = open(template_fasta_file, 'w')
    template_fasta_fp.write('>%s\n' % left_template_name)
    template_fasta_fp.write('%s\n' % left_template_seq)
    template_fasta_fp.write('>%s\n' % right_template_name)
    template_fasta_fp.write('%s\n' % right_template_seq)
    template_fasta_fp.close()

    first_round_paf_file = os.path.join(out_dir, 'first_round.paf')
    preset = tk.get_preset_for_minimap2(platform)
    cmd = f'{minimap2} -A 2 -B 5 -c --eqx -t {num_threads} -x {preset} {template_fasta_file} {in_fastq_file} > {first_round_paf_file} 2> /dev/null'

    tk.eprint('NOTICE: first round alignment')
    tk.run_system_cmd(cmd)
    tk.eprint('NOTICE: first round alignment finished')

    first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set = first_round_estimation_from_paf(first_round_paf_file, len(left_anchor_seq), len(right_anchor_seq), repeat_unit)
    
    return first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set

def first_round_estimation_from_paf(first_round_paf_file, left_anchor_len, right_anchor_len, repeat_unit):

    first_round_repeat_count_dict = dict()
    potential_repeat_region_dict = dict()
    bad_reads_set = set()
    first_round_paf_fp = open(first_round_paf_file, 'r')

    read_paf_list = list()
    while 1:
        line = first_round_paf_fp.readline()
        if not line: break
        line = line.strip()
        if not line: continue

        col_list = line.strip().split('\t')
        paf = tk.PAF(col_list)

        if len(read_paf_list) == 0 or paf.qname == read_paf_list[0].qname:
            read_paf_list.append(paf)
        else:
            first_round_estimation_for1read(read_paf_list, left_anchor_len, right_anchor_len, repeat_unit, first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set)
            read_paf_list.clear()
            read_paf_list.append(paf)

    first_round_paf_fp.close()
    
    first_round_estimation_for1read(read_paf_list, left_anchor_len, right_anchor_len, repeat_unit, first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set)

    return first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set

def first_round_estimation_for1read(read_paf_list, left_anchor_len, right_anchor_len, repeat_unit, first_round_repeat_count_dict, potential_repeat_region_dict, bad_reads_set):

    repeat_unit_size = len(repeat_unit)
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
            sys.exit()
    
    if len(left_anchor_paf_list) != 1 and len(right_anchor_paf_list) != 1 : return
    
    outer_end1 = -1
    outer_end2 = -1
    if len(left_anchor_paf_list) == 1:
        # left anchor is good
        left_paf  = left_anchor_paf_list[0]
        left_paf_num_repeats = int((left_paf.tend - left_boundary_pos)/repeat_unit_size) + 1
        repeat_size1 = calculate_repeat_size_from_exact_match(left_paf.cigar, left_paf.tstart, left_boundary_pos, repeat_unit_size)
        outer_end1 = left_paf.qend
        if left_paf.strand == '-': 
            outer_end1 = left_paf.qlen - outer_end1
        
    if len(right_anchor_paf_list) == 1:
        # right anchor is good
        right_paf = right_anchor_paf_list[0]
        right_paf_num_repeats = int((right_paf.tend - right_boundary_pos)/repeat_unit_size) + 1
        repeat_size2 = calculate_repeat_size_from_exact_match(right_paf.cigar, right_paf.tstart, right_boundary_pos, repeat_unit_size)
        outer_end2 = right_paf.qend
        if right_paf.strand == '-': 
            outer_end2 = right_paf.qlen - outer_end2

    
    if len(left_anchor_paf_list) != 1:
        # left anchor is bad
        first_round_repeat_count_dict[right_paf.qname] = (repeat_size2, right_paf_num_repeats)
        return
    
    if len(right_anchor_paf_list) != 1:
        # rigth anchor is bad
        first_round_repeat_count_dict[left_paf.qname] = (repeat_size1, left_paf_num_repeats)
        return

    if left_paf.strand == right_paf.strand: return

    if (left_paf.strand == '+' and outer_end1 < outer_end2) or (left_paf.strand == '-' and outer_end1 > outer_end2):
        bad_reads_set.add(left_paf.qname)
        return
    
    max_num_repeat_units = max(left_paf_num_repeats, right_paf_num_repeats)+5
    min_num_repeat_units = min(repeat_size1, repeat_size2)
    if left_paf.strand == '-': 
        outer_end1, outer_end2 = tk.switch_two_numbers(outer_end1, outer_end2)

    flank_len = 100
    a = outer_end2 - flank_len
    if a < 0: a = 0
    b = outer_end1 + flank_len
    if b > left_paf.qlen: b = left_paf.qlen
    potential_repeat_region_dict[left_paf.qname] = (a, b)

    if left_paf.qname not in first_round_repeat_count_dict:
        first_round_repeat_count_dict[left_paf.qname] = (min_num_repeat_units, max_num_repeat_units)
    else:
        tk.eprint(f'ERROR! read name already exists in first_round_repeat_count_dict: {left_paf.qname}')
        sys.exit()

    return 

def extract_anchor_sequence(repeat_chrom_seq, repeat_start, repeat_end, max_anchor_len):
    
    anchor_start = repeat_start - max_anchor_len
    anchor_end = repeat_start
    if anchor_start < 0: anchor_start = 0
    left_anchor_seq = repeat_chrom_seq[anchor_start:anchor_end]

    anchor_start = repeat_end
    anchor_end = anchor_start + max_anchor_len
    if anchor_end > len(repeat_chrom_seq):
        anchor_end = len(repeat_chrom_seq)
    right_anchor_seq = repeat_chrom_seq[anchor_start:anchor_end]
    
    return left_anchor_seq, right_anchor_seq

def calculate_repeat_size_from_exact_match(cigar, tstart, ref_repeat_start_pos, repeat_unit_size):
    repeat_size = 0
    cigar_opr_list, cigar_opr_len_list = tk.analysis_cigar_string(cigar)
    current_ref_pos = tstart
    for i in range(0, len(cigar_opr_list)):
        cigar_opr = cigar_opr_list[i]
        cigar_opr_len = cigar_opr_len_list[i]
        if cigar_opr == '=':  # match
            if current_ref_pos >= ref_repeat_start_pos:
                repeat_size += int(cigar_opr_len/repeat_unit_size)
            elif current_ref_pos + cigar_opr_len >= ref_repeat_start_pos:
                overlap_len = current_ref_pos + cigar_opr_len - ref_repeat_start_pos
                if overlap_len > 0: repeat_size += int(overlap_len/repeat_unit_size)
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'X': # mismatch
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'I': # insertion
            pass
        elif cigar_opr == 'D': # deletion
            current_ref_pos   += cigar_opr_len
        else:
            tk.eprint('ERROR! unsupported CIGAR operation: %s' % cigar_opr)
            sys.exit()

    return repeat_size

def analysis_repeat_region(repeat_region):

    a = repeat_region.split(':')
    if len(a) != 2:
        tk.eprint('ERROR! --repeat_region should be in the format of chr:start-end (e.g. chr4:3074876-3074939)!')
        sys.exit(1)
    
    repeat_chr = a[0]

    b = a[1].split('-')
    if len(b) != 2:
        tk.eprint('ERROR! --repeat_region should be in the format of chr:start-end (e.g. chr4:3074876-3074939)!')
        sys.exit(1)

    repeat_start = int(b[0])
    repeat_end = int(b[1])

    return repeat_chr, repeat_start, repeat_end

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
        return '%d\t%d\t%d\t%d\t%d' % (self.allele_id, self.num_reads, self.predicted_repeat_size, self.min_repeat_size, self.max_repeat_size)


class FastqChunk:
    def __init__(self, min_size, max_size, fn, fp):
        self.min_size = int(min_size)
        self.max_size = int(max_size)
        self.fn = fn
        self.fp = fp


if __name__ == '__main__':
    main()
