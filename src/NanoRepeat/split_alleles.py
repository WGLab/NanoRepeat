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
import random
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import math
import gzip

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
matplotlib.rcParams['font.family'] = "sans-serif"

from NanoRepeat import tk
from NanoRepeat.repeat_region import *

class Allele:
    def __init__(self):
        self.gmm_mean1 = None
        self.gmm_mean2 = None
        self.gmm_sd1 = None
        self.gmm_sd2 = None
        self.readname_list = []
        self.repeat1_size_list = []
        self.repeat2_size_list = []
        self.repeat1_median_size = None
        self.repeat2_median_size = None
        self.probability_list = []
        self.num_reads = None
        self.allele_frequency = None
        self.confidence_list = []

        self.gmm_min1 = None
        self.gmm_min2 = None
        self.gmm_max1 = None
        self.gmm_max2 = None

class Readinfo:
    def __init__(self, readname):

        self.readname = readname
        self.label = -1
        self.repeat_size1 = -1
        self.repeat_size2 = -1
        self.confidence = 1

def simulate_reads(read_repeat_count_list, error_rate):
    simulated_read_repeat_count_list = read_repeat_count_list * 100
    for i in range(0, len(simulated_read_repeat_count_list)):
        std = error_rate * (10 + simulated_read_repeat_count_list[i])
        random_error = random.gauss(0, std)
        simulated_read_repeat_count_list[i] += random_error
    return simulated_read_repeat_count_list

def interval_has_overlap(interval1, interval2):
    start = min(interval1[1], interval2[1])
    end = max(interval1[0], interval2[0])
    if end - start <= 0: 
        return True
    else:
        return False

def get_outlier_cutoff_from_list(repeat_count_list):

    if len(repeat_count_list) == 0:
        tk.eprint('ERROR! len(repeat_count_list) == 0')
        tk.eprint('Please report this bug!')
        sys.exit(1)

    mean = np.mean(repeat_count_list)
    std = np.std(repeat_count_list)

    min_repeat_count_cutoff = mean - 3 * std
    if min_repeat_count_cutoff < 0: min_repeat_count_cutoff = 0
    max_repeat_count_cutoff = mean + 3 * std

    return min_repeat_count_cutoff, max_repeat_count_cutoff

def analysis_outlier_1d(read_repeat_count_dict):

    repeat_count_list = list()
    for readname in read_repeat_count_dict:
        repeat_count = read_repeat_count_dict[readname]
        repeat_count_list.append(repeat_count)
    min_repeat_count_cutoff, max_repeat_count_cutoff = get_outlier_cutoff_from_list(repeat_count_list)
 
    return min_repeat_count_cutoff, max_repeat_count_cutoff

def remove_outlier_reads_2d(read_repeat_joint_count_dict):

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

    return readname_list, read_repeat_count_list

def remove_outlier_reads_1d(read_repeat_count_dict):

    min_count, max_count  = analysis_outlier_1d(read_repeat_count_dict)

    readname_list = list() # readnames after removing outliers
    read_repeat_count_list = list()

    for readname in read_repeat_count_dict:
        repeat_count = read_repeat_count_dict[readname]
        if repeat_count < min_count or repeat_count > max_count: continue 
        readname_list.append(readname)
        read_repeat_count_list.append(repeat_count)

    return readname_list, read_repeat_count_list

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

def auto_GMM_1d(X, max_num_components, cov_type, max_mutual_overlap):

    for n in range(2, max_num_components + 1):
        gmm = GaussianMixture(n_components=n, covariance_type=cov_type, n_init=10).fit(X)
        # gmm.means_.shape = (n, 1) (n_components, n_variable)
        # gmm.covariances_.shape = (n, 1) (n_components, n_variable)
        for i in range(0, n):
            for j in range(i+1, n):
                component_i_mean = gmm.means_[i]
                component_j_mean = gmm.means_[j]
                component_i_cov  = gmm.covariances_[i]
                component_j_cov  = gmm.covariances_[j]
                
                component_i_sd = max(1.0, math.sqrt(component_i_cov))
                component_j_sd = max(1.0, math.sqrt(component_j_cov))

                a1 = 1.0-max_mutual_overlap
                a2 = max_mutual_overlap
                interval_i = (norm.isf(a1, component_i_mean, component_i_sd), norm.isf(a2, component_i_mean, component_i_sd))
                interval_j = (norm.isf(a1, component_j_mean, component_j_sd), norm.isf(a2, component_j_mean, component_j_sd))

                if interval_has_overlap(interval_i, interval_j):
                    best_num_components = n - 1
                    gmm = GaussianMixture(n_components=best_num_components, covariance_type=cov_type, n_init=10).fit(X)
                    return best_num_components, gmm

    
    best_num_components = max_num_components
    gmm = GaussianMixture(n_components=best_num_components, covariance_type=cov_type, n_init=10).fit(X)
    return best_num_components, gmm


def auto_GMM_2d(X, max_num_components, cov_type, max_mutual_overlap):
    for n in range(1, max_num_components+1):
        gmm = GaussianMixture(n_components=n, covariance_type=cov_type, n_init=10).fit(X)
        # gmm.means_.shape = (n, 2) (n_components, n_variable)
        # gmm.covariances_.shape = (n, 2) (n_components, n_variable)
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
                if interval_has_overlap(interval_i1, interval_j1) and interval_has_overlap(interval_i2, interval_j2):
                    best_num_components = n - 1
                    gmm = GaussianMixture(n_components=best_num_components, covariance_type=cov_type, n_init=10).fit(X)
                    return best_num_components, gmm

    
    best_num_components = max_num_components
    gmm = GaussianMixture(n_components=best_num_components, covariance_type=cov_type, n_init=10).fit(X)
    return best_num_components, gmm

def create_allele_list_1d(best_n_components, final_gmm, readname_list, read_repeat_count_array, read_repeat_count_dict, probability_cutoff):

    read_label_list = list(final_gmm.predict(read_repeat_count_array))
    proba2darray = final_gmm.predict_proba(read_repeat_count_array)

    allele_list = []
    for i in range(0, best_n_components):
        allele = Allele()
        allele.gmm_mean1 = final_gmm.means_[i]
        allele.gmm_sd1 = math.sqrt(final_gmm.covariances_[i])
        allele_list.append(allele)

    assert len(readname_list) == len(read_label_list)
    
    for i in range(0, len(readname_list)):
        readname = readname_list[i]
        repeat_size1 = read_repeat_count_dict[readname]
        assert repeat_size1 == read_repeat_count_array[i][0]
        read_label = read_label_list[i]
        probability = proba2darray[i][read_label]
        allele_list[read_label].readname_list.append(readname)
        allele_list[read_label].repeat1_size_list.append(repeat_size1)
        allele_list[read_label].probability_list.append(probability)

    for read_label in range(0, len(allele_list)):
        allele_list[read_label].num_reads = len(allele_list[read_label].readname_list)
        if allele_list[read_label].num_reads == 0:
            allele_list[read_label].repeat1_median_size = 0
            allele_list[read_label].gmm_min1 = 0
            allele_list[read_label].gmm_max1
            continue

        allele_list[read_label].repeat1_median_size = int(np.median(allele_list[read_label].repeat1_size_list)+0.5)
        allele_list[read_label].gmm_min1 = allele_list[read_label].gmm_mean1 - 2 * allele_list[read_label].gmm_sd1
        allele_list[read_label].gmm_max1 = allele_list[read_label].gmm_mean1 + 2 * allele_list[read_label].gmm_sd1

    for read_label in range(0, len(allele_list)):
        allele_list[read_label].confidence_list = []
        for i in range(0, len(allele_list[read_label].readname_list)):
            confidence = 'HIGH'
            if allele_list[read_label].probability_list[i] < probability_cutoff: 
                confidence = 'LOW'
            if allele_list[read_label].repeat1_size_list[i] < allele_list[read_label].gmm_min1 or allele_list[read_label].repeat1_size_list[i] > allele_list[read_label].gmm_max1:
                confidence = 'LOW'
            allele_list[read_label].confidence_list.append(confidence)

    allele_list.sort(key = lambda allele:allele.num_reads)

    while allele_list[0].num_reads == 0:
        allele_list.pop(0)
    
    return allele_list

def create_allele_list_2d(best_n_components, final_gmm, readname_list, read_repeat_count_array, read_repeat_joint_count_dict, probability_cutoff):
    
    read_label_list = list(final_gmm.predict(read_repeat_count_array))
    proba2darray = final_gmm.predict_proba(read_repeat_count_array)

    allele_list = []
    for i in range(0, best_n_components):
        allele = Allele()
        allele.gmm_mean1 = final_gmm.means_[i][0]
        allele.gmm_mean2 = final_gmm.means_[i][1]
        allele.gmm_sd1 = math.sqrt(final_gmm.covariances_[i][0])
        allele.gmm_sd2 = math.sqrt(final_gmm.covariances_[i][1])
        allele_list.append(allele)

    assert len(readname_list) == len(read_label_list)
    for i in range(0, len(readname_list)):
        readname = readname_list[i]
        repeat_size1, repeat_size2 = read_repeat_joint_count_dict[readname]
        read_label = read_label_list[i]
        probability = proba2darray[i][read_label]
        allele_list[read_label].readname_list.append(readname)
        allele_list[read_label].repeat1_size_list.append(repeat_size1)
        allele_list[read_label].repeat2_size_list.append(repeat_size2)
        allele_list[read_label].probability_list.append(probability)
    
    for read_label in range(0, len(allele_list)):
        allele_list[read_label].num_reads = len(allele_list[read_label].readname_list)
        if allele_list[read_label].num_reads == 0:
            allele_list[read_label].repeat1_median_size = 0
            allele_list[read_label].repeat2_median_size = 0
            allele_list[read_label].gmm_min1 = 0
            allele_list[read_label].gmm_min2 = 0
            allele_list[read_label].gmm_max1 = 0
            allele_list[read_label].gmm_max2 = 0
            continue
        
        allele_list[read_label].repeat1_median_size = int(np.median(allele_list[read_label].repeat1_size_list)+0.5)
        allele_list[read_label].repeat2_median_size = int(np.median(allele_list[read_label].repeat2_size_list)+0.5)
        allele_list[read_label].gmm_min1 = allele_list[read_label].gmm_mean1 - 2 * allele_list[read_label].gmm_sd1
        allele_list[read_label].gmm_min2 = allele_list[read_label].gmm_mean2 - 2 * allele_list[read_label].gmm_sd2
        allele_list[read_label].gmm_max1 = allele_list[read_label].gmm_mean1 + 2 * allele_list[read_label].gmm_sd1
        allele_list[read_label].gmm_max2 = allele_list[read_label].gmm_mean2 + 2 * allele_list[read_label].gmm_sd2

    for read_label in range(0, len(allele_list)):
        allele_list[read_label].confidence_list = []
        for i in range(0, len(allele_list[read_label].readname_list)):
            confidence = 'HIGH'
            if allele_list[read_label].probability_list[i] < probability_cutoff: 
                confidence = 'LOW'
            if allele_list[read_label].repeat1_size_list[i] < allele_list[read_label].gmm_min1 or allele_list[read_label].repeat1_size_list[i] > allele_list[read_label].gmm_max1:
                confidence = 'LOW'
            if allele_list[read_label].repeat2_size_list[i] < allele_list[read_label].gmm_min2 or allele_list[read_label].repeat2_size_list[i] > allele_list[read_label].gmm_max2:
                confidence = 'LOW'
            allele_list[read_label].confidence_list.append(confidence)

    allele_list.sort(key = lambda allele:allele.num_reads)

    while allele_list[0].num_reads == 0:
        allele_list.pop(0)

    return allele_list



def create_readinfo_dict_from_allele_list(allele_list, dimension):

    if dimension not in [1, 2]:
        tk.eprint('ERROR! dimension must be 1 or 2! Please report this bug.')
        sys.exit(1)
    
    readinfo_dict = dict()
    for read_label in range(0, len(allele_list)):
        for i in range(0, len(allele_list[read_label].readname_list)):
            readname = allele_list[read_label].readname_list[i]
            readinfo = Readinfo(readname)
            readinfo.label = read_label
            readinfo.repeat_size1 = allele_list[read_label].repeat1_size_list[i]
            if dimension == 2:
                readinfo.repeat_size2 = allele_list[read_label].repeat2_size_list[i]

            readinfo.confidence = allele_list[read_label].confidence_list[i]
            readinfo_dict[readname] = readinfo

    return readinfo_dict

def output_phasing_results_1d(repeat_region:RepeatRegion, allele_list):

    if repeat_region.no_details == False:
        out_phasing_file = repeat_region.out_prefix + '.phased_reads.txt'
        out_phasing_f = open(out_phasing_file, 'w')
        header  = f'##RepeatRegion={repeat_region.to_unique_id()}\n'
        header += f'#Read_Name\tAllele_ID\tPhasing_Confidence\tRepeat_Size\n'
        out_phasing_f.write(header)
        out = ''

    for label in range(0,len(allele_list)):
        allele_id = label + 1
        allele = allele_list[label]

        for i in range(0, len(allele.readname_list)):
            readname = allele.readname_list[i]
            repeat_size1 = allele.repeat1_size_list[i]
            confidence = allele.confidence_list[i]

            if repeat_region.no_details == False:
                out += f'{readname}\t{allele_id}\t{confidence}\t{repeat_size1:.1f}\n'

            if readname not in repeat_region.results.quantified_read_dict:
                repeat_region.results.quantified_read_dict[readname] = QuantifiedRead()
                repeat_region.results.quantified_read_dict[readname].read_name = readname

            repeat_region.results.quantified_read_dict[readname].repeat_size1 = repeat_size1
            repeat_region.results.quantified_read_dict[readname].allele_id = allele_id
            repeat_region.results.quantified_read_dict[readname].phasing_confidence = confidence

    if repeat_region.no_details == False:
        out_phasing_f.write(out)
        out_phasing_f.close()

    return

def output_phasing_results_2d(allele_list, repeat1_id, repeat2_id, in_fastq_file, out_prefix):

    out_phasing_file = out_prefix + '.phased_reads.txt'
    out_phasing_f = open(out_phasing_file, 'w')
    header  = f'##Input_FASTQ={in_fastq_file}\n'
    header += f'#Read_Name\tAllele_ID\tPhasing_Confidence\t{repeat1_id}.Repeat_Size\t{repeat2_id}.Repeat_Size\n'
    out_phasing_f.write(header)
    out = ''
    for label in range(0,len(allele_list)):
        allele_id = label + 1
        allele = allele_list[label]
        for i in range(0, len(allele.readname_list)):
            readname = allele.readname_list[i]
            repeat_size1 = allele.repeat1_size_list[i]
            repeat_size2 = allele.repeat2_size_list[i]
            confidence = allele.confidence_list[i]
            out += f'{readname}\t{allele_id}\t{confidence}\t{repeat_size1:.1f}\t{repeat_size2:.1f}\n'
    
    out_phasing_f.write(out)
    out_phasing_f.close()

    return


def output_phased_fastq(in_fastq_file, readinfo_dict, num_alleles, out_prefix):

    out_allele_fastq_file_list = list()
    for label in range(0, num_alleles):
        allele_id = label + 1
        out_allele_fastq_file = out_prefix + '.allele%d.fastq' % (allele_id)
        out_allele_fastq_file_list.append(out_allele_fastq_file)

    out_allele_fastq_f_list = list()
    for i in range(0, len(out_allele_fastq_file_list)):
        out_allele_fastq_f = open(out_allele_fastq_file_list[i], 'w')
        out_allele_fastq_f_list.append(out_allele_fastq_f)
    
    if '.gz' == in_fastq_file[-3:]:
        in_fastq_f = gzip.open(in_fastq_file, 'rt')
    else:
        in_fastq_f = open(in_fastq_file, 'rt')

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
        if readname not in readinfo_dict: continue
        if readinfo_dict[readname].confidence != 'HIGH': continue
    
        label = readinfo_dict[readname].label
        out_allele_fastq_f_list[label].write(line1 + line2 + line3 + line4)

    in_fastq_f.close()

    for i in range(0, len(out_allele_fastq_f_list)):
        out_allele_fastq_f_list[i].close()

    return

def output_summary_file_1d(repeat_region:RepeatRegion, allele_list, num_removed_reads, out_prefix):

    num_alleles = len(allele_list)
    if repeat_region.no_details == False:
        
        out_summray_file = out_prefix + '.summary.txt'
        filebasename = os.path.split(out_summray_file)[1]
        out_summray_f = open(out_summray_file, 'w')
        summary_info  = f'Summary_file={filebasename}\tRepeat_Region={repeat_region.to_unique_id()}'
        summary_info += f'\tMethod=GMM'
        summary_info += f'\tNum_Alleles={num_alleles}'
        summary_info += f'\tNum_Removed_Reads={num_removed_reads}'

    for label in range(0, num_alleles):
        allele_id = label + 1
        if repeat_region.no_details == False:
            summary_info += f'\tAllele{allele_id}_Num_Reads={allele_list[label].num_reads}'
            summary_info += f'\tAllele{allele_id}_Repeat_Size={allele_list[label].repeat1_median_size}'

        quantified_allele = QuantifiedAllele()
        quantified_allele.repeat_size1 = allele_list[label].repeat1_median_size
        quantified_allele.num_supp_reads = allele_list[label].num_reads
        repeat_region.results.quantified_allele_list.append(quantified_allele)

        
    if repeat_region.no_details == False:
        summary_info += '\n'
        out_summray_f.write(summary_info)
        out_summray_f.close()

    return

def output_summary_file_2d(in_fastq_file, allele_list, repeat1_id, repeat2_id, num_removed_reads, out_prefix):

    num_alleles = len(allele_list)
    out_summray_file = out_prefix + '.summary.txt'
    out_summray_f = open(out_summray_file, 'w')
    summary_info  = f'Input_FASTQ\t{in_fastq_file}\n'
    summary_info += f'Method\t2D-GMM\n'
    summary_info += f'Num_Alleles\t{num_alleles}\n'
    summary_info += f'Num_Removed_Reads\t{num_removed_reads}\n'

    for label in range(0, num_alleles):
        allele_id = label + 1
        summary_info += f'Allele{allele_id}_Num_Reads\t{allele_list[label].num_reads}\n'
        summary_info += f'Allele{allele_id}_{repeat1_id}.Repeat_Size\t{allele_list[label].repeat1_median_size}\n'
        summary_info += f'Allele{allele_id}_{repeat2_id}.Repeat_Size\t{allele_list[label].repeat2_median_size}\n'

    out_summray_f.write(summary_info)
    out_summray_f.close()

    return

def output_repeat_size_1d(repeat_region:RepeatRegion):

    if repeat_region.no_details == False:
        repeat_size_file =  f'{repeat_region.out_prefix}.repeat_size.txt'
        repeat_size_f = open(repeat_size_file, 'w')
        header  = f'##Repeat_Region={repeat_region.to_unique_id()}\n'
        header += f'#Read_Name\tRepeat_Size\n'
        repeat_size_f.write(header)
    
    for read_name in repeat_region.read_dict:
        repeat_size = repeat_region.read_dict[read_name].round3_repeat_size
        if repeat_size != None:
            if repeat_region.no_details == False:
                repeat_size_f.write(f'{read_name}\t{repeat_size:.1f}\n')
            if read_name not in repeat_region.results.quantified_read_dict:
                repeat_region.results.quantified_read_dict[read_name] = QuantifiedRead()
                repeat_region.results.quantified_read_dict[read_name].read_name = read_name
                repeat_region.results.quantified_read_dict[read_name].repeat_size1 = repeat_size

    if repeat_region.no_details == False:
        repeat_size_f.close()

    return

def output_repeat_size_2d (in_fastq_file, repeat1_id, repeat2_id, out_prefix, repeat1_count_dict, repeat2_count_dict):

    repeat_size_file = f'{out_prefix}.repeat_size.txt'
    repeat_size_f = open(repeat_size_file, 'w')
    header = f'##Input_FASTQ={in_fastq_file}\n'
    header += f'#Read_Name\t{repeat1_id}.Repeat_Size\t{repeat2_id}.Repeat_Size\n'
    repeat_size_f.write(header)

    all_readname_list = []
    for readname in repeat1_count_dict:
        all_readname_list.append(readname)
    for readname in repeat2_count_dict:
        all_readname_list.append(readname)

    all_readname_set = set(all_readname_list)
    read_repeat_joint_count_dict = dict()
    read_repeat_joint_count_list = []
    for readname in all_readname_set:
        if readname in repeat1_count_dict: 
            size1 = repeat1_count_dict[readname]
        else:
            size1 = 'N.A.'

        if readname in repeat2_count_dict:
            size2 = repeat2_count_dict[readname]
        else:
            size2 = 'N.A.'
        
        read_repeat_joint_count_dict[readname] = (size1, size2)
        tup = (readname, size1, size2)
        read_repeat_joint_count_list.append(tup)
    
    read_repeat_joint_count_list.sort(key = lambda x:x[1])
    for tup in read_repeat_joint_count_list:
        readname, size1, size2 = tup
        repeat_size_f.write(f'{readname}\t{size1:.1f}\t{size2:.1f}\n')
    repeat_size_f.close()

    tk.eprint(f'NOTICE: Repeat size file is here: {repeat_size_file}')
    return read_repeat_joint_count_dict


def plot_repeat_counts_1d(readinfo_dict, allele_list, repeat_id, out_prefix):

    hist_figure_file = out_prefix + '.hist.png'
    num_alleles = len(allele_list)
    x_list = list()

    x_2d_list = [0] * num_alleles
    for i in range(0, num_alleles):
        x_2d_list[i] = list()

    for readname in readinfo_dict:
        readinfo = readinfo_dict[readname]

        x = readinfo.repeat_size1

        x_list.append(x)

        x_2d_list[readinfo.label].append(x)

    xmin = int(min(x_list))
    xmax = int(max(x_list))+1

    if xmax - xmin <= 200:
        b1 = range(xmin - 1, xmax + 2)
    else:
        b1 = range(xmin - 1, xmax + 2, int(float(xmax - xmin)/200.0 + 0.5))

    repeat1_predicted_size_list = []
    for label in range(0, len(allele_list)):
        repeat1_predicted_size_list.append(allele_list[label].repeat1_median_size)

    plot_hist1d(x_2d_list, b1, repeat_id, repeat1_predicted_size_list, allele_list, hist_figure_file)

    return

def plot_repeat_counts_2d(readinfo_dict, allele_list, repeat1_id, repeat2_id, out_prefix):

    hist2d_figure_file = out_prefix + '.hist2d.png'
    repeat1_hist_figure_file  = out_prefix + '.%s.hist.png' % (repeat1_id)
    repeat2_hist_figure_file  = out_prefix + '.%s.hist.png' % (repeat2_id)

    num_alleles = len(allele_list)
    x_list = list()
    y_list = list()

    x_2d_list = [0] * num_alleles
    y_2d_list = [0] * num_alleles
    for i in range(0, num_alleles):
        x_2d_list[i] = list()
        y_2d_list[i] = list()

    for readname in readinfo_dict:
        readinfo = readinfo_dict[readname]

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
    
    repeat1_predicted_size_list = []
    repeat2_predicted_size_list = []
    for label in range(0, len(allele_list)):
        repeat1_predicted_size_list.append(allele_list[label].repeat1_median_size)
        repeat2_predicted_size_list.append(allele_list[label].repeat2_median_size)

    plot_hist2d(x_list, y_list, b1, b2, repeat1_id, repeat2_id, allele_list, hist2d_figure_file)
    plot_hist1d(x_2d_list, b1, repeat1_id, repeat1_predicted_size_list, allele_list, repeat1_hist_figure_file)
    plot_hist1d(y_2d_list, b2, repeat2_id, repeat2_predicted_size_list, allele_list, repeat2_hist_figure_file)

    return 

def get_1d_max_min_from_allele_list(allele_list):

    gmm_min_list = []
    gmm_max_list = []

    for allele in allele_list:
        gmm_min_list.append(allele.gmm_min1)
        gmm_max_list.append(allele.gmm_max1)

    xmin = min(gmm_min_list)
    xmax = max(gmm_max_list)

    xmax = int(((xmax * 1)/10.0 +1 )* 10)
    xmin = int(((xmin)/10.0 -1 )* 10)

    if xmin < 10: xmin = 0

    return xmin, xmax

def get_2d_max_min_from_allele_list(allele_list):

    assert(len(allele_list) > 0)

    gmm_min1_list = []
    gmm_max1_list = []
    gmm_min2_list = []
    gmm_max2_list = []

    for allele in allele_list:
        gmm_min1_list.append(allele.gmm_min1)
        gmm_min2_list.append(allele.gmm_min2)
        gmm_max1_list.append(allele.gmm_max1)
        gmm_max2_list.append(allele.gmm_max2)

    xmin = min(gmm_min1_list)
    xmax = max(gmm_max1_list)
    ymin = min(gmm_min2_list)
    ymax = max(gmm_max2_list)

    xmax = int(((xmax * 1)/10.0 +1 )* 10)
    ymax = int(((ymax * 1)/10.0 +1 )* 10)

    xmin = int(((xmin)/10.0 -1 )* 10)
    ymin = int(((ymin)/10.0 -1 )* 10)

    if xmin < 10: xmin = 0
    if ymin < 10: ymin = 0
    
    return xmin, xmax, ymin, ymax


def plot_hist1d(x_2d_list, b, repeat_id, predicted_size_list, allele_list, out_file):

    plt.figure (figsize=(6, 4))
    
    max_x = 0
    for x in x_2d_list:
        if len(x) == 0: continue
        plt.hist(x, bins = b)
        if max_x < max(x):
            max_x = max(x)

    for repeat_size in predicted_size_list:
        plt.axvline(x=repeat_size+0.5, color = 'grey', linestyle = ':')

    if len(repeat_id) > 30:
        repeat_id = repeat_id[0:30] + '...'
    
    plt.title('Repeat size distribution (%s)' % repeat_id)
    plt.xlabel('repeat size')
    plt.ylabel('number of reads')
    
    xmin, xmax = get_1d_max_min_from_allele_list(allele_list)
    
    plt.xlim(xmin, xmax)

    plt.savefig(out_file, dpi=300)
    plt.close('all')

    return

def plot_hist2d(x_list, y_list, b1, b2, repeat_id1, repeat_id2, allele_list, out_file):

    plt.figure (figsize=(8, 4))
    plt.hist2d (x_list, y_list, [b1, b2], cmap = 'binary')
    plt.colorbar()
    plt.title('2D histogram of repeat size')
    plt.xlabel('repeat size (%s)' % repeat_id1)
    plt.ylabel('repeat size (%s)' % repeat_id2)

    xmin, xmax, ymin, ymax = get_2d_max_min_from_allele_list(allele_list)
   
    '''
    xmin = min(x_list)
    xmax = max(x_list)
    ymin = min(y_list)
    ymax = max(y_list)
    '''
    
    if xmax - xmin < 20: 
        delta = (20 - (xmax - xmin))/2
        xmin -= delta
        if xmin < 0: xmin = 0
        xmax = xmin + 20
    if ymax - ymin < 20: 
        delta = (20-(ymax-ymin))/2
        ymin -= delta
        if ymin < 0: ymin = 0
        ymax = ymin + 20
    
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.savefig(out_file, dpi=300)
    plt.close('all')

    return

def scatter_plot_with_contour_2d (read_repeat_joint_count_dict, final_gmm, score_cut_off, repeat1_id, repeat2_id, allele_list, out_prefix):

    scatter_plot_file = out_prefix + '.scatter.png'

    xlabel = '%s repeat size' % (repeat1_id)
    ylabel = '%s repeat size' % (repeat2_id)

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

    figure_min_x, figure_max_x, figure_min_y, figure_max_y = get_2d_max_min_from_allele_list(allele_list)

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