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
import string
import sys
import shutil
import numpy as np
from multiprocessing import Process
from typing import List

from NanoRepeat import tk
from NanoRepeat.repeat_region import *
from NanoRepeat.paf import *
from NanoRepeat.split_alleles import *
from NanoRepeat.repeat_region import *

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

def refine_repeat_region_in_ref(minimap2:string, repeat_region:RepeatRegion):

    mid_ref_seq            = repeat_region.mid_ref_seq
    repeat_unit_seq        = repeat_region.repeat_unit_seq
    repeat_region_location = repeat_region.to_invertal()
    temp_out_dir           = repeat_region.temp_out_dir
    out_prefix             = repeat_region.out_prefix
    error_message_file     = out_prefix + '.error_message.log'
    mid_ref_seq_len        = len(mid_ref_seq)

    buffer_size = 10
    if mid_ref_seq_len < buffer_size: return

    num_units = int(float(mid_ref_seq_len) / len(repeat_unit_seq))+1
    pure_repeat_seq = repeat_unit_seq * num_units

    ref_check_paf_file   = os.path.join(temp_out_dir, 'ref_aligned_to_pure_repeats.paf')
    mid_ref_seq_file     = os.path.join(temp_out_dir, 'mid_ref_seq.fasta')
    pure_repeat_seq_file = os.path.join(temp_out_dir, 'pure_repeat_seq.fasta')

    write_seq_to_fasta(f'mid_ref_seq_{repeat_region_location}', mid_ref_seq, mid_ref_seq_file)
    write_seq_to_fasta('pure_repeat_seq', pure_repeat_seq, pure_repeat_seq_file)

    score_parameters = '-x ava-ont -z30 -k3 -w2 -m1 -n2 -s10'

    cmd = f'{minimap2} {score_parameters} -f 0 --cs  {pure_repeat_seq_file} {mid_ref_seq_file} > {ref_check_paf_file}'
    tk.run_system_cmd(cmd)

    ref_check_paf_f = open(ref_check_paf_file, 'r')
    lines = list(ref_check_paf_f)
    ref_check_paf_f.close()

    paf_list = []
    for line in lines:
        col_list = line.strip().split('\t')
        paf = PAF(col_list)
        paf_list.append(paf)
    
    if len(paf_list) == 0:
        repeat_region.ref_has_issue = True
        
        error_message = f'ERROR! The given repeat motif ({repeat_region.repeat_unit_seq}) was not found in the reference sequence ({repeat_region_location}). Please check your input!\n'
        error_message += f'The repeat motif is: {repeat_unit_seq}\n'
        error_message += f'The sequence inside {repeat_region_location} is: {mid_ref_seq}\n'
        
        tk.eprint(f'current repeat: {repeat_region.to_unique_id()}\n')
        tk.eprint(error_message)
        error_message_f = open(error_message_file, 'w')
        error_message_f.write(error_message)
        error_message_f.close()

        return

    max_align_score     = paf_list[0].align_score
    max_align_score_idx = 0
    for i in range(1, len(paf_list)):
        if paf_list[i].align_score > max_align_score:
            max_align_score = paf_list[i].align_score
            max_align_score_idx = i
    
    best_paf  = paf_list[max_align_score_idx]
    align_len = best_paf.qend - best_paf.qstart

    if best_paf.align_score < 0:
        repeat_region.ref_has_issue = True
        error_message  = f'ERROR! In the reference sequence, the interval ({repeat_region_location}) has too many mismatches with the repeat motif. Please check your input.\n\n'
        error_message += f'The repeat motif is: {repeat_unit_seq}\n'
        error_message += f'The sequence inside {repeat_region_location} is:\n'
        error_message += mid_ref_seq[0:best_paf.qstart] + '---' + mid_ref_seq[best_paf.qstart:best_paf.qend] + '---' + mid_ref_seq[best_paf.qend:] + '\n'
        
        tk.eprint(f'current repeat: {repeat_region.to_unique_id()}\n')
        tk.eprint(error_message)
        error_message_f = open(error_message_file, 'w')
        error_message_f.write(error_message)
        error_message_f.close()
        return
        
    if mid_ref_seq_len - align_len > buffer_size and float(align_len)/mid_ref_seq_len < 0.8:
        repeat_region.ref_has_issue = True
        
        corrected_start_pos = repeat_region.start_pos + best_paf.qstart
        corrected_end_pos   = corrected_start_pos + align_len
        
        error_message = f'ERROR! In the reference sequence, the interval ({repeat_region_location}) included non-repeat sequences or repeats that are not the given motif. Please check your input.\n'
        error_message += f'The repeat motif is: {repeat_unit_seq}\n'
        error_message += f'The sequence inside {repeat_region_location} is:\n'
        error_message += mid_ref_seq[0:best_paf.qstart] + '---' + mid_ref_seq[best_paf.qstart:best_paf.qend] + '---' + mid_ref_seq[best_paf.qend:] + '\n'
        error_message += f'\n\nWe suggest you change the chromosome positions (in the `--repeat_region_bed` file) to:\n{repeat_region.chrom}\t{corrected_start_pos}\t{corrected_end_pos}\t{repeat_unit_seq}\n'
        
        tk.eprint(f'current repeat: {repeat_region.to_unique_id()}\n')
        tk.eprint(error_message)
        error_message_f = open(error_message_file, 'w')
        error_message_f.write(error_message)
        error_message_f.close()

    return

def write_seq_to_fasta(seq_name, seq, out_file):

    out_f = open(out_file, 'w')
    out_f.write(f'>{seq_name}\n')
    out_f.write(f'{seq}\n')
    out_f.close()

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

    if len(read_paf_list) == 0: return

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

def find_anchor_locations_in_reads(minimap2:string, data_type:string, repeat_region:RepeatRegion, num_cpu:int):

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

    preset = tk.get_preset_for_minimap2(data_type)
    cmd = f'{minimap2} -c -t {num_cpu} {preset} {template_fasta_file} {repeat_region.region_fq_file} > {anchor_locations_paf_file}'
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

        repeat_region.read_core_seq_dict[read_name] = core_sequence
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


def round1_and_round2_estimation(minimap2:string, data_type:string, repeat_region: RepeatRegion, num_cpu: int):

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
    
    preset = tk.get_preset_for_minimap2(data_type)
    cmd = f'{minimap2} -c -t {num_cpu} {preset} -f 0.0 {round1_fasta_file} {repeat_region.core_seq_fq_file} > {round1_paf_file}'
    tk.run_system_cmd(cmd)
    round1_paf_f = open(round1_paf_file, 'r')
    lines = list(round1_paf_f)
    round1_paf_f.close()


    round2_repeat_size_dict = dict()
    
    for line in lines:
        col_list = line.strip().split('\t')
        paf = PAF(col_list)

        if paf.tstart <= len(repeat_region.left_anchor_seq) and paf.tend >= len(repeat_region.left_anchor_seq):
            read_name = paf.qname
            round2_repeat_size = float(paf.tend - len(repeat_region.left_anchor_seq))/len(repeat_region.repeat_unit_seq)
            align_score = paf.align_score
            if read_name not in round2_repeat_size_dict:
                round2_repeat_size_dict[read_name] = []
            round2_repeat_size_dict[read_name].append((align_score, round2_repeat_size))
            
    for read_name in round2_repeat_size_dict:
        round2_repeat_size_dict[read_name].sort(key = lambda x:x[0], reverse=True)
        round2_repeat_size = round2_repeat_size_dict[read_name][0][1]
        repeat_region.read_dict[read_name].round2_repeat_size = round2_repeat_size

    round2_repeat_size_file = os.path.join(repeat_region.temp_out_dir, 'round2_repeat_size.txt')
    round2_repeat_size_f = open(round2_repeat_size_file, 'w')
    for read_name in repeat_region.read_dict:
        round2_repeat_size_f.write(f'{read_name}\t{repeat_region.read_dict[read_name].round2_repeat_size}\n')
    
    round2_repeat_size_f.close()
    
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

def round3_estimation(minimap2:string, data_type:string, fast_mode, repeat_region:RepeatRegion, num_cpu:int):
    
    round3_paf_file = os.path.join(repeat_region.temp_out_dir, 'round3.paf')
    round3_paf_f = open(round3_paf_file, 'w')
    round3_paf_f.close()
    
    procs = []

    for proc_id in range(0, num_cpu):
        proc = Process(target=round3_align_1process, args=(proc_id, num_cpu, fast_mode, repeat_region, data_type, minimap2))
        procs.append(proc)
        proc.start()
        
    for proc_id in range(0, num_cpu):
        procs[proc_id].join()
        proc_out_paf_file = os.path.join(repeat_region.temp_out_dir, f'round3_align.process{proc_id}.paf')
        cmd = f'cat {proc_out_paf_file} >> {round3_paf_file}'
        tk.run_system_cmd(cmd)
        os.remove(proc_out_paf_file)
    
    round3_estimation_from_paf(repeat_region, round3_paf_file)

def round3_align_1process(proc_id:int, num_cpu:int, fast_mode, repeat_region:RepeatRegion, data_type:string, minimap2:string):
    
    read_id = -1
    out_paf_file = os.path.join(repeat_region.temp_out_dir, f'round3_align.process{proc_id}.paf')
    out_paf_f = open(out_paf_file, 'w')
    out_paf_f.close()
    
    for read_name in repeat_region.read_dict:
        read = repeat_region.read_dict[read_name]
        round2_repeat_size = read.round2_repeat_size
        if round2_repeat_size == None: continue
        read_id += 1
        if read_id % num_cpu != proc_id: continue
        
        buffer = max(15, int(round2_repeat_size * 0.05))
        if buffer > 150: buffer = 150
        
        if fast_mode:
            buffer = 15
        
        max_template_repeat_size = int(round2_repeat_size + buffer)
        min_template_repeat_size = int(round2_repeat_size - buffer)
        
        if min_template_repeat_size < 0: min_template_repeat_size = 0
        
        template_fasta_file  = os.path.join(repeat_region.temp_out_dir, f'round3_reference.read{read_id}.fasta')
        template_fasta_fp    = open(template_fasta_file, 'w')
        for repeat_size in range(min_template_repeat_size, max_template_repeat_size+1):
            template_seq = repeat_region.left_anchor_seq + repeat_region.repeat_unit_seq * repeat_size + repeat_region.right_anchor_seq
            template_fasta_fp.write('>%s\n' % repeat_size)
            template_fasta_fp.write('%s\n' % template_seq)

        template_fasta_fp.close()
        
        read_seq = repeat_region.read_core_seq_dict[read_name].strip()
        read_fasta_file = os.path.join(repeat_region.temp_out_dir, f'round3_input.read{read_id}.fasta')
        
        read_fasta_f = open(read_fasta_file, 'w')
        read_fasta_f.write(f'>{read_name}\n')
        read_fasta_f.write(read_seq + '\n')
        read_fasta_f.close()

        preset = tk.get_preset_for_minimap2(data_type)
        cmd = f'{minimap2} {preset} -f 0.0 -N 100 -c --eqx -t 1 {template_fasta_file} {read_fasta_file} >> {out_paf_file}'
        tk.run_system_cmd(cmd)

        os.remove(template_fasta_file)
        os.remove(read_fasta_file)

    return

def remove_noisy_reads_1d(allele_list, ploidy):

    tk.eprint('NOTICE: Trying to remove noisy reads')
    allele_list.sort(key = lambda allele:allele.num_reads)
    num_removed_reads = 0
    while len(allele_list) > ploidy and len(allele_list) >= 2:
        if allele_list[0].num_reads * 1.5 <= allele_list[-ploidy].num_reads:
            num_removed_reads += allele_list[0].num_reads 
            allele_list.pop(0)
        else:
            break

    return allele_list, num_removed_reads


def split_allele_using_gmm_1d(repeat_region:RepeatRegion, ploidy, error_rate, max_mutual_overlap, max_num_components, remove_noisy_reads):

    in_fastq_file = repeat_region.region_fq_file
    out_prefix = repeat_region.out_prefix

    if ploidy < 1:
        tk.eprint('ERROR! ploidy must be >= 1 !')
        sys.exit()


    read_repeat_count_dict = dict()
    for read_name in repeat_region.read_dict:
        if repeat_region.read_dict[read_name].round3_repeat_size != None:
            read_repeat_count_dict[read_name] = repeat_region.read_dict[read_name].round3_repeat_size

    if len(read_repeat_count_dict) == 0:
        tk.eprint(f'ERROR! No reads were found for repeat region: {repeat_region.to_outfile_prefix()}')
        return

    if len(read_repeat_count_dict) == 1:
        tk.eprint(f'ERROR! No enough reads for phasing. Repeat region is: {repeat_region.to_outfile_prefix()}')
        return

    probability_cutoff = 0.95
    cov_type = 'diag'
    dimension = 1
    
    readname_list, read_repeat_count_list = remove_outlier_reads_1d(read_repeat_count_dict)
    read_repeat_count_array = np.array(read_repeat_count_list).reshape(-1, 1)

    simulated_read_repeat_count_list = simulate_reads(read_repeat_count_list, error_rate)
    simualted_read_repeat_count_array = np.array(simulated_read_repeat_count_list).reshape(-1, dimension)

    best_n_components, final_gmm = auto_GMM_1d (simualted_read_repeat_count_array, max_num_components, cov_type, max_mutual_overlap)
    tk.eprint(f'NOTCIE: Number of alleles={best_n_components}')
    
    allele_list = create_allele_list_1d(best_n_components, final_gmm, readname_list, read_repeat_count_array, read_repeat_count_dict, probability_cutoff)
    
    num_removed_reads = 0
    if remove_noisy_reads == True and len(allele_list) > ploidy:
        allele_list, num_removed_reads = remove_noisy_reads_1d(allele_list, ploidy)
        tk.eprint(f'NOTICE: There are {len(allele_list)} alleles after removing noisy reads')
    
    best_n_components = len(allele_list)
    allele_list.sort(key = lambda allele:allele.gmm_mean1)

    readinfo_dict = create_readinfo_dict_from_allele_list(allele_list, dimension)

    tk.eprint('NOTICE: Writing phasing results...')
    output_phasing_results_1d(repeat_region, allele_list)
    
    output_phased_fastq(in_fastq_file, readinfo_dict, best_n_components, out_prefix)

    tk.eprint('NOTICE: Writing summary file...')
    repeat_region.results.num_alleles = len(allele_list)
    output_summary_file_1d(repeat_region, allele_list, num_removed_reads, out_prefix)

    tk.eprint('NOTICE: Plotting figures...')
    if repeat_region.no_details == False:
        plot_repeat_counts_1d(readinfo_dict, allele_list, repeat_region.to_unique_id(), out_prefix)

    return
    
def quantify1repeat_from_bam(input_args, error_rate, in_bam_file:string, ref_fasta_dict:dict, repeat_region:RepeatRegion):

    def _clean_and_exit(input_args, repeat_region:RepeatRegion):

        if input_args.save_temp_files == False:
            shutil.rmtree(repeat_region.temp_out_dir)

        if repeat_region.no_details == True:
            shutil.rmtree(f'{input_args.out_prefix}.details')

        return repeat_region.results

    temp_out_dir = f'{input_args.out_prefix}.NanoRepeat_temp_dir.{repeat_region.to_outfile_prefix()}'
    out_dir = f'{input_args.out_prefix}.details/{repeat_region.chrom}'
    os.makedirs(temp_out_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    repeat_region.temp_out_dir = temp_out_dir
    repeat_region.out_prefix = f'{out_dir}/{repeat_region.to_outfile_prefix()}'
    repeat_region.anchor_len = input_args.anchor_len

    # extract reads from bam file
    repeat_region.region_fq_file = os.path.join(temp_out_dir, f'{repeat_region.to_outfile_prefix()}.fastq')
    region_bam_file = f'{repeat_region.out_prefix}.sorted.bam'
    cmd = f'{input_args.samtools} view --reference {input_args.ref_fasta} -hb {in_bam_file} {repeat_region.to_invertal(flank_dist=repeat_region.anchor_len)} > {region_bam_file}'
    tk.run_system_cmd(cmd)
    
    cmd = f'{input_args.samtools} index {region_bam_file}'
    tk.run_system_cmd(cmd)
    
    cmd = f'{input_args.samtools} fastq {region_bam_file} > {repeat_region.region_fq_file}'
    tk.run_system_cmd(cmd)

    fastq_file_size = os.path.getsize(repeat_region.region_fq_file)
    
    if fastq_file_size == 0:
        tk.eprint(f'WARNING! No reads were found in repeat region: {repeat_region.to_outfile_prefix()}')
        return _clean_and_exit(input_args, repeat_region)

    # extract ref sequence
    extract_ref_sequence(ref_fasta_dict, repeat_region)
    
    refine_repeat_region_in_ref(input_args.minimap2, repeat_region)
    
    if repeat_region.ref_has_issue == True and input_args.save_temp_files == False:
        return _clean_and_exit(input_args, repeat_region)
    
    tk.eprint('NOTICE: Step 1: finding anchor location in reads')
    find_anchor_locations_in_reads(input_args.minimap2, input_args.data_type, repeat_region, input_args.num_cpu)
    tk.eprint('NOTICE: Step 1 finished')

    # make core sequence fastq
    make_core_seq_fastq(repeat_region)

    tk.eprint('NOTICE: Step 2: round 1 and round 2 estimation')
    round1_and_round2_estimation(input_args.minimap2, input_args.data_type, repeat_region, input_args.num_cpu)
    tk.eprint('NOTICE: Step 2 finished')

    tk.eprint('NOTICE: Step 3: round 3 estimation')
    round3_estimation(input_args.minimap2, input_args.data_type, input_args.fast_mode, repeat_region, input_args.num_cpu)
    tk.eprint('NOTICE: Step 3 finished')

    tk.eprint('NOTICE: Writing to repeat size file...')
    output_repeat_size_1d(repeat_region)

    tk.eprint('NOTICE: Step 4: phasing reads using GMM')
    split_allele_using_gmm_1d(repeat_region, input_args.ploidy, error_rate, input_args.max_mutual_overlap, input_args.max_num_components, input_args.remove_noisy_reads)

    return _clean_and_exit(input_args, repeat_region)

def nanoRepeat_bam (input_args, in_bam_file:string):
    
    # ont, ont_sup, ont_q20, clr, hifi

    if input_args.data_type == 'ont' or 'clr':
        error_rate = 0.07
    elif input_args.data_type == 'ont_sup':
        error_rate = 0.04
    elif input_args.data_type == 'ont_q20':
        error_rate = 0.03
    elif input_args.data_type == 'hifi':
        error_rate = 0.02
    else:
        tk.eprint(f'ERROR! unknown data type: {input_args.data_type}')
        sys.exit(1)
    
    tk.eprint(f'NOTICE: Reading repeat region file: {input_args.repeat_region_bed}')
    repeat_region_list = read_repeat_region_file(input_args.repeat_region_bed, input_args.no_details)

    tk.eprint(f'NOTICE: Reading reference fasta file: {input_args.ref_fasta}')
    ref_fasta_dict = tk.fasta_file2dict(input_args.ref_fasta)
    
    out_tsv_file = f'{input_args.out_prefix}.NanoRepeat_output.tsv'
    out_tsv_f = open(out_tsv_file, 'w')
    for repeat_region in repeat_region_list:
        tk.eprint(f'NOTICE: Quantifying repeat: {repeat_region.to_outfile_prefix()}')
        results = quantify1repeat_from_bam(input_args, error_rate, in_bam_file, ref_fasta_dict, repeat_region)
        out_tsv_f.write(f'{repeat_region.to_tab_invertal()}\t{repeat_region.repeat_unit_seq}\t')
        num_alleles = len(results.quantified_allele_list)
        out_tsv_f.write(f'{num_alleles}\t{results.max_repeat_size1()}\t{results.min_repeat_size1()}\t{results.allele_summary()}\t{results.read_summary()}\n')
        
    out_tsv_f.close()
    tk.eprint('NOTICE: Program finished.')

    return

