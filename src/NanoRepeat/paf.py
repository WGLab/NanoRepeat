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
from NanoRepeat import tk

class PAF:
    def __init__(self, col_list):
        if len(col_list) < 12:
            tk.eprint('ERROR! number of columns should be >= 12 in a PAF file!')
            tk.eprint('The corrupted line is: %s' % '\t'.join(col_list))
            sys.exit()
        
        self.qname, self.qlen, self.qstart, self.qend = col_list[0:4]
        self.strand = col_list[4]
        self.tname, self.tlen, self.tstart, self.tend = col_list[5:9]
        self.n_match, self.align_len, self.mapq = col_list[9:12]

        self.qlen      = int(self.qlen)
        self.qstart    = int(self.qstart)
        self.qend      = int(self.qend)
        self.tlen      = int(self.tlen)
        self.tstart    = int(self.tstart)
        self.tend      = int(self.tend)
        self.n_match   = int(self.n_match)
        self.align_len = int(self.align_len)
        self.mapq      = int(self.mapq)
        self.is_primary= False
        self.align_score = -1
        self.cigar = ''
        for col in col_list[12:]:
            if col[0:5] == 'AS:i:': 
                self.align_score = int(col[5:])
            elif col[0:5] == 'cg:Z:':
                self.cigar = col[5:]
            elif col == 'tp:A:P':
                self.is_primary = True
            elif col == 'tp:A:S':
                self.is_primary = False

        if self.strand != '+' and self.strand != '-':
            tk.eprint('ERROR! unknown strand: %s' % self.strand)
            sys.exit()
        
        if self.strand == '-':
            new_qstart  = self.qlen - self.qend
            new_qend    = self.qlen - self.qstart
            self.qstart = new_qstart
            self.qend   = new_qend
        
        self.tscore = 0

    def output_core(self):
        return f'{self.qname}\t{self.strand}\t{self.qlen}\t{self.qstart}\t{self.qend}\t{self.tname}\t{self.tlen}\t{self.tstart}\t{self.tend}'
        

        