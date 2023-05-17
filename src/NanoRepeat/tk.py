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

import gzip
from datetime import datetime
import distutils
from distutils import spawn
from timeit import repeat

TimeFormat = '%m/%d/%Y %H:%M:%S'

class ReadError:
    def __init__(self, num_match = 0, num_mismatch = 1000000, num_ins = 1000000, num_del = 1000000, align_score = 0):
        self.num_match = num_match
        self.num_mismatch = num_mismatch
        self.num_ins = num_ins
        self.num_del = num_del
        self.score = align_score
    
    def num_edit_bases(self):
        return self.num_mismatch + self.num_ins + self.num_del

### IO ###

def read_list_file(in_list_file, abspath = False):
    out_list = list()
    in_list_fp = open(in_list_file, 'r')
    lines = list(in_list_fp)
    in_list_fp.close()

    for line in lines:
        line = line.strip()
        if abspath == True:
            line = os.path.abspath(line)
            out_list.append(line)
        else:
            out_list.append(line)
    return out_list

def gzopen(in_file, mode = 'rt'):

    if '.gz' == in_file[-3:] or '.bgz' == in_file[-4:]:
        in_fp = gzip.open(in_file, 'rt')
    else:
        in_fp = open(in_file, 'rt')

    return in_fp

def create_dir(dir):

    os.makedirs(dir, exist_ok=True)

    return

def check_input_file_exists(in_file):

    if os.path.exists(in_file) == False:
        eprint('ERROR! the file was not found: %s' % in_file)
        sys.exit()

    return

def check_output_file_exists(out_file):

    if os.path.exists(out_file) == False or os.path.getsize(out_file) == 0:
        return False
    else:
        return True

def eprint(message):
    sys.stderr.write('[' + datetime.now().strftime(TimeFormat) + '] ' + message + '\n')
    return

def get_file_prefix(input_file):

    return os.path.splitext(os.path.split(input_file)[1])[0]

### FASTQ/FASTA ###
    
def count_fastq (in_fastq_file):

    n_reads = 0

    in_fastq_fp = gzopen(in_fastq_file)
    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        n_reads += 1

    in_fastq_fp.close()

    return n_reads

def fasta_file2dict(fasta_file):

    fasta_fp = gzopen(fasta_file)
    fasta_dict = dict()

    curr_name = ''
    curr_seq = ''

    while 1:
        line = fasta_fp.readline()
        if not line: break
        line = line.strip()
        if not line: continue 
        if line[0] ==  '>':
            if len(curr_seq) > 0 and len(curr_name) > 0:
                assert curr_name not in fasta_dict
                fasta_dict[curr_name] = curr_seq
            curr_name = line[1:].split()[0]
            curr_seq = ''
            continue

        curr_seq += line
     
    fasta_fp.close()

    if len(curr_seq) > 0 and len(curr_name) > 0:
        fasta_dict[curr_name] = curr_seq
    
    return fasta_dict

def read_fasta_file(fasta_file):

    fasta_fp = gzopen(fasta_file)

    fasta_name_list = list()
    fasta_seq_list = list()
    curr_name = ''
    curr_seq = ''

    while 1:
        line = fasta_fp.readline()
        if not line: break
        line = line.strip()
        if not line: continue 
        if line[0] ==  '>':
            if len(curr_seq) > 0 and len(curr_name) > 0:
                fasta_name_list.append(curr_name)
                fasta_seq_list.append(curr_seq)
            curr_name = line[1:].split()[0]
            curr_seq = ''
            continue

        curr_seq += line
     
    fasta_fp.close()

    if len(curr_seq) > 0 and len(curr_name) > 0:
        fasta_name_list.append(curr_name)
        fasta_seq_list.append(curr_seq)
    
    return fasta_name_list, fasta_seq_list


def read_one_chr_from_fasta_file(fasta_file, target_chr):

    if '.gz' == fasta_file[-3:]:
        fasta_fp = gzip.open(fasta_file, 'rt')
    else:
        fasta_fp = open(fasta_file, 'rt')

    target_seq = ''

    status = 's' # skip
    while 1:
        line = fasta_fp.readline()
        if not line: break
        line = line.strip()
        if not line: continue
        if line[0] ==  '>':
            new_contig_name = line[1:].split()[0]
            if new_contig_name == target_chr:
                status = 'r'
                continue
            else:
                if len(target_seq) > 0:
                    status = 'd'
                    break
                else:
                    status = 's'
                    continue

        if status == 'd': 
            break
        elif status == 's':
            continue
        elif status == 'r':
            target_seq += line
     
    fasta_fp.close()
    
    return target_seq


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

def split_fastq(in_fastq_file_list, num_out_file, out_prefix):

    out_readname_set = set() # used to skip duplicated reads
    out_fastq_file_list = list()
    out_fastq_fp_list = list()
    for i in range(0, num_out_file):
        out_fastq_file = out_prefix + '.part%d.fastq' % i
        out_fastq_file_list.append(out_fastq_file)
        out_fastq_fp = open(out_fastq_file, 'w')
        out_fastq_fp_list.append(out_fastq_fp)

    for in_fastq_file in in_fastq_file_list:
        in_fastq_fp = gzopen(in_fastq_file)

        i = 0
        status = 0
        while 1:
            if status == 0:
                line1 = in_fastq_fp.readline()
            line2 = in_fastq_fp.readline()
            line3 = in_fastq_fp.readline()
            line4 = in_fastq_fp.readline()
            if not line1: break
            if not line2: break
            if not line3: break
            if not line4: break

            if line1[0] != '@' or len(line2) != len(line4) or line3.strip() != '+':
                eprint('ERROR! Bad fastq file: %s' % in_fastq_file)
                break
                

            i += 1
            file_id = i % num_out_file
            readname = line1.strip().split()[0][1:]
            if readname not in out_readname_set:
                out_fastq_fp_list[file_id].write(line1)
                out_fastq_fp_list[file_id].write(line2)
                out_fastq_fp_list[file_id].write(line3)
                out_fastq_fp_list[file_id].write(line4)

            out_readname_set.add(readname)

        in_fastq_fp.close()

    for out_fastq_fp in out_fastq_fp_list:
        out_fastq_fp.close()

    return out_fastq_file_list

def extract_fastq_tail_seq(in_fastq_file, read_tail_length, left_tail_fastq_file, right_tail_fastq_file):

    in_fastq_fp = gzopen(in_fastq_file)
    fq_left_tail_fp  = open(left_tail_fastq_file, 'w')
    fq_right_tail_fp = open(right_tail_fastq_file, 'w')

    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        if line1[0] != '@' or len(line2) != len(line4) or line3.strip() != '+':
            eprint('ERROR! bad fastq record: ')
            eprint(line1 + line2 + line3 + line4)
            sys.exit()
        
        read_seq  = line2.strip()
        read_qual = line4.strip()

        left_tail_seq  = read_seq[0:read_tail_length]
        left_tail_qual = read_qual[0:read_tail_length]

        right_tail_seq  = rev_comp(read_seq[-read_tail_length:])
        right_tail_qual = ''.join(reversed(read_qual[-read_tail_length:]))

        fq_left_tail_fp.write(line1)
        fq_left_tail_fp.write(left_tail_seq + '\n')
        fq_left_tail_fp.write(line3)
        fq_left_tail_fp.write(left_tail_qual + '\n')

        fq_right_tail_fp.write(line1)
        fq_right_tail_fp.write(right_tail_seq + '\n')
        fq_right_tail_fp.write(line3)
        fq_right_tail_fp.write(right_tail_qual + '\n')

    in_fastq_fp.close()
    fq_left_tail_fp.close()
    fq_right_tail_fp.close()

    return


def rev_comp(seq):

    complement  = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rev_seq = ''.join(reversed(seq))

    rev_comp_seq = ''
    for b in rev_seq:
        rev_comp_seq += complement[b]

    return rev_comp_seq



def format_parameters(input_para):

    input_para = input_para.strip('`')
    input_para = input_para.strip('\'')
    input_para = input_para.strip('\"')
    input_para = input_para.lower()

    return input_para

def compute_overlap_len(start1, end1, start2, end2):
    max_start = max(start1, start2)
    min_end   = min(end1, end2)

    overlap_len = min_end - max_start
    return max(0, overlap_len)



### SAM/BAM/Alignment ###


def analysis_cigar_string (cigar):
    
    opr_set   = set(['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P', 'M'])
    digit_set = set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])

    cigar_opr_list = list()
    cigar_opr_len_list = list()
    
    length_str = ''
    for i in range(0, len(cigar)):
        c = cigar[i]
        if c in digit_set:
            length_str += c
        elif c in opr_set:
            cigar_opr_list.append(c)
            cigar_opr_len_list.append(int(length_str))
            length_str = ''
        else:
            eprint('ERROR: unknown CIGAR operation: %s' % c)
            sys.exit()
    
    return cigar_opr_list, cigar_opr_len_list



def calculate_repeat_size_from_exact_match(cigar, tstart, ref_repeat_start_pos, repeat_unit_size):
    repeat_size = 0
    cigar_opr_list, cigar_opr_len_list = analysis_cigar_string(cigar)
    current_ref_pos = tstart

    for i in range(0, len(cigar_opr_list)):
        cigar_opr = cigar_opr_list[i]
        cigar_opr_len = cigar_opr_len_list[i]
        if cigar_opr == '=':  # match
            if current_ref_pos >= ref_repeat_start_pos:
                repeat_size += int(cigar_opr_len/repeat_unit_size)
            elif current_ref_pos + cigar_opr_len >= ref_repeat_start_pos:
                overlap_len = current_ref_pos + cigar_opr_len - ref_repeat_start_pos
                if overlap_len > 0: 
                    repeat_size += int(overlap_len/repeat_unit_size)
        
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'X': # mismatch
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'I': # insertion
            pass
        elif cigar_opr == 'D': # deletion
            current_ref_pos   += cigar_opr_len
        else:
            eprint('ERROR! unsupported CIGAR operation: %s' % cigar_opr)
            sys.exit()

    return repeat_size


def target_region_alignment_stats_from_cigar(cigar, tstart, tend, ref_region_start, ref_region_end):

    # tstart is the aligned start position in the paf file
    # ref_region_start and ref_region_end are the target region in the reference

    if len(cigar) == 0:
        eprint('ERROR! cigar string is empty!')
        sys.exit()
    
    matching_score = 2
    mismatching_penalty = -4
    gap_open_penalty = -4
    gap_ext_penalty = -2

    num_match = 0
    num_mismatch = 0
    num_ins = 0
    num_del = 0
    align_score = 0
    cigar_opr_list, cigar_opr_len_list = analysis_cigar_string(cigar)
    if len(cigar_opr_list) != len(cigar_opr_len_list):
        eprint('ERROR: len(cigar_opr_list) != len(cigar_opr_len_list) This could be a bug.')
        sys.exit()
        
    current_ref_pos = tstart
  
    for i in range(0, len(cigar_opr_list)):
        cigar_opr = cigar_opr_list[i]
        cigar_opr_len = cigar_opr_len_list[i]
        if cigar_opr == '=':  # match
            overlap_len = compute_overlap_len(current_ref_pos, current_ref_pos+cigar_opr_len, ref_region_start, ref_region_end)
            if overlap_len > 0: 
                num_match += overlap_len
                align_score += overlap_len * matching_score
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'X': # mismatch
            overlap_len = compute_overlap_len(current_ref_pos, current_ref_pos+cigar_opr_len, ref_region_start, ref_region_end)
            if overlap_len > 0: 
                num_mismatch += overlap_len
                align_score += overlap_len * mismatching_penalty
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'I': # insertion
            if current_ref_pos > ref_region_start and current_ref_pos < ref_region_end - 1:
                num_ins += cigar_opr_len
                align_score += gap_open_penalty + (cigar_opr_len-1) * gap_ext_penalty
        elif cigar_opr == 'D': # deletion
            overlap_len = compute_overlap_len(current_ref_pos, current_ref_pos+cigar_opr_len, ref_region_start, ref_region_end)
            if overlap_len > 0:
                num_del += overlap_len
                align_score += gap_open_penalty + (overlap_len-1) * gap_ext_penalty
            current_ref_pos += cigar_opr_len
        elif cigar_opr == 'S':
            continue
        else:
            eprint('ERROR: unsupported cigar operation: %s' % cigar_opr)
            sys.exit()
        if current_ref_pos > ref_region_end: break

    if tend < ref_region_end:
        num_mismatch += ref_region_end - tend
    
    if tstart > ref_region_start:
        num_mismatch += tstart - ref_region_start

    read_error = ReadError(num_match, num_mismatch, num_ins, num_del, align_score)
    return read_error

def get_preset_for_minimap2(data_type):
    if data_type == 'ont':
        preset = ' -x map-ont '
    elif data_type == 'ont_sup':
        preset = ' -x map-ont '
    elif data_type == 'ont_q20':
        preset = ' -x map-ont '
    elif data_type == 'clr':
        preset = ' -x map-ont '
    elif data_type == 'hifi':
        preset = ' -x map-ont '
    else:
        eprint(f'ERROR: Unknown data type: {data_type}\n')
        sys.exit(1)
    
    return preset

minimap2_mid_para   = ' -k 5 -w 3 -n 1 -m 10 -s 40 '
minimap2_short_para = ' -k 3 -w 2 -n 1 -m 10 -s 40 '
### system ###

def switch_two_numbers(a, b):
    return b, a

def switch_two_objects(a, b):
    return b, a

def run_system_cmd(cmd):
    eprint(f'NOTICE: Running command: {cmd}')  
    ret = os.system(cmd)
    if ret != 0: 
        eprint('ERROR: Failed to run command: %s' % cmd)
        eprint('Return value is: %d' % ret)
        sys.exit()
    return

def find_executable_path(name):
    path = distutils.spawn.find_executable(name)
    return path

def rm(path):
    if os.path.exists(path):
        try:
            os.remove(path)
        except:
            eprint(f'WARNING: Failed to remove file: {path}')
    return
