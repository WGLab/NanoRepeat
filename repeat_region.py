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

## objects ##
class Read:
    def __init__(self, readname:str, readseq:str):
        self.read_name = readname
        self.read_seq = readseq
        self.core_seq = None
        self.is_good_for_analysis = None
        self.init_repeat_size = None
    

class RepeatRegion:
    def __init__(self, line = None):
        self.left_anchor_seq = None
        self.right_anchor_seq = None
        
        self.left_anchor_len = None
        self.right_anchor_len = None

        self.anchor_len = None
        self.mid_ref_seq = None

        self.region_fq_file = None
        self.region_fasta_file = None # template

        self.chrom = None
        self.start_pos = None
        self.end_pos = None
        self.repeat_unit_seq = None
        self.max_repeat_size = None

        self.temp_out_dir = None
        self.final_out_dir = None

        self.temp_file_list = []
    
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

    def to_unique_id(self):
        return f'{self.chrom}-{self.start_pos}-{self.end_pos}-{self.repeat_unit_seq}'

def read_repeat_region_file(repeat_region_file):
    repeat_region_list = []
    repeat_region_f = open(repeat_region_file, 'r')
    lines = list(repeat_region_f)
    repeat_region_f.close()
    for line in lines:
        repeat_region = RepeatRegion(line)
        repeat_region_list.append(repeat_region)

    return repeat_region_list
