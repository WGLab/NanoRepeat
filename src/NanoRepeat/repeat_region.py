#!/usr/bin/env python3

'''
Copyright (c) 2020-2023 Children's Hospital of Philadelphias
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

## objects ##
class Read:
    def __init__(self):
        self.read_name = None
        self.full_read_len = None
        self.core_seq = None
        self.init_repeat_size = None
        self.left_anchor_is_good = False
        self.right_anchor_is_good = False
        self.both_anchors_are_good = False
        self.core_seq_start_pos = None
        self.core_seq_end_pos = None
        self.mid_seq_start_pos = None
        self.mid_seq_end_pos = None
        self.left_anchor_paf = None
        self.right_anchor_paf = None
        self.dist_between_anchors = None
        self.seq_between_anchors = None
        self.strand = None
        self.left_buffer_len = None
        self.right_buffer_len = None
        self.round1_repeat_size = None
        self.round2_repeat_size = None
        self.round3_repeat_size = None

class QuantifiedAllele:
    def __init__(self):
        self.num_supp_reads = '*'
        self.repeat_size1 = '*'
        self.repeat_size2 = '*'

class QuantifiedRead:
    def __init__(self):
        self.read_name = '*'
        self.repeat_size1 = -1
        self.repeat_size2 = -1
        self.allele_id = -1
        self.phasing_confidence = -1

class Result:
    def __init__(self):
        self.quantified_allele_list = []
        self.quantified_read_dict = dict()

    def allele_summary(self):
        summary = 'Allele_Repeat_Size;Allele_Num_Support_Reads'
        for i in range(0, len(self.quantified_allele_list)):
            quantified_allele = self.quantified_allele_list[i]
            summary += f'|{quantified_allele.repeat_size1};{quantified_allele.num_supp_reads}'

        return summary
    
    def read_summary(self):
        summary = 'Read_Name;Read_Repeat_Size;Read_Allele_ID;PhasingConfidence'
        for read_name in self.quantified_read_dict:
            quantified_read = self.quantified_read_dict[read_name]
            summary += f'|{quantified_read.read_name};{quantified_read.repeat_size1};{quantified_read.allele_id};{quantified_read.phasing_confidence}'

        return summary
    
    def max_repeat_size1(self):
        if len(self.quantified_allele_list) == 0:
            return -1
        
        max_size = self.quantified_allele_list[0].repeat_size1

        for i in range(1, len(self.quantified_allele_list)):
            if self.quantified_allele_list[i].repeat_size1 > max_size:
                max_size = self.quantified_allele_list[i].repeat_size1
        return max_size
    
    def min_repeat_size1(self):
        if len(self.quantified_allele_list) == 0:
            return -1
        
        min_size = self.quantified_allele_list[0].repeat_size1

        for i in range(1, len(self.quantified_allele_list)):
            if self.quantified_allele_list[i].repeat_size1 < min_size:
                min_size = self.quantified_allele_list[i].repeat_size1
        return min_size
    

class RepeatRegion:
    def __init__(self, line = None, no_details=False):
        self.left_anchor_seq = None
        self.right_anchor_seq = None
        
        self.left_anchor_len = None
        self.right_anchor_len = None

        self.anchor_len = None
        self.mid_ref_seq = None

        self.region_fq_file = None
        self.region_fasta_file = None # template
        self.core_seq_fq_file = None
        self.mid_seq_fq_file = None
        self.chrom = None
        self.start_pos = None
        self.end_pos = None
        self.repeat_unit_seq = None
        self.max_repeat_size = None

        self.temp_out_dir = None
        self.final_out_dir = None
        self.out_prefix = None

        self.temp_file_list = []
        self.ref_has_issue = False
        self.read_dict = dict()
        self.read_core_seq_dict = dict()
        self.buffer_len = None
        self.brute_force_repeat_count_dict = dict()

        self.no_details = no_details
        self.results = Result()

        if line != None:
            col_list = line.strip().split('\t')
            if len(col_list) < 4:
                sys.stderr.write('ERROR! the repeat region bed file should be tab-delimited and have 4 columns: chrom, start_position, end_position, repeat_unit. start_position and end_position should be 0-based')
                sys.exit(1)
            
            self.chrom, self.start_pos, self.end_pos, self.repeat_unit_seq = col_list[0:4]
            self.start_pos = int(self.start_pos)
            self.end_pos = int(self.end_pos)
    
    def to_invertal(self, flank_dist = 0):
        assert (flank_dist >= 0)
        start_pos = self.start_pos - flank_dist
        end_pos = self.end_pos + flank_dist
        if start_pos < 0: start_pos = 0
        return f'{self.chrom}:{start_pos}-{end_pos}'

    def to_tab_invertal(self, flank_dist = 0):
        assert (flank_dist >= 0)
        start_pos = self.start_pos - flank_dist
        end_pos = self.end_pos + flank_dist
        if start_pos < 0: start_pos = 0
        return f'{self.chrom}\t{start_pos}\t{end_pos}'
    
    def to_unique_id(self):
        return f'{self.chrom}-{self.start_pos}-{self.end_pos}-{self.repeat_unit_seq}'

    def to_outfile_prefix(self):
        if len(self.repeat_unit_seq) < 30:
            seq = self.repeat_unit_seq
        else:
            seq = self.repeat_unit_seq[0:20] + '....' + self.repeat_unit_seq[-6:]

        return f'{self.chrom}-{self.start_pos}-{self.end_pos}-{seq}'
    
def read_repeat_region_file(repeat_region_file, no_details):
    repeat_region_list = []
    repeat_region_f = open(repeat_region_file, 'r')
    lines = list(repeat_region_f)
    repeat_region_f.close()
    for line in lines:
        repeat_region = RepeatRegion(line, no_details)
        repeat_region_list.append(repeat_region)

    return repeat_region_list
