#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import functools
import mimetypes
import gzip
from functools import partial
from Bio.SeqIO.QualityIO import FastqGeneralIterator
"""
Filter a Illumina FASTQ file based on index sequence.

Reads a Illumina FASTQ file and compares the sequence index in the
`sample number` position of the sequence identifier to a supplied sequence
index. Entries that match the sequence index are filtered into the *filtered
file* (if any) and entries that don't match are filtered into the *unfiltered
file* (if any). Displays the count of total, filtered and unfiltered reads.
Matching with mismatches (`-m` parameter), and gzip compression for input
(detected on the basis of file extension) and output (specified using `-c`
parameter) are supported.
"""
_PROGRAM_NAME = 'filter_illumina_index'
# -------------------------------------------------------------------------------
# Author:       Tet Woo Lee
#
# Created:      2018-12-13
# Copyright:    © 2018 Tet Woo Lee
# Licence:      GPLv3
#
# Dependencies: Biopython, tested on v1.72
# -------------------------------------------------------------------------------

_PROGRAM_VERSION = '1.0.1'
# -------------------------------------------------------------------------------
# ### Change log
#
# version 1.0.1 2018-12-19
# : Speed up number of mismatches calculation
#
# version 1.0 2018-12-14
# : Minor updates for PyPi and conda packaging
#
# version 1.0.dev1 2018-12-13
# : First working version
# -------------------------------------------------------------------------------


# INITIALISATION

parser = argparse.ArgumentParser(prog=_PROGRAM_NAME,
                                 description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_required_named = parser.add_argument_group('required named arguments')
parser.add_argument('--version', action='version',
                    version='{} {}'.format(_PROGRAM_NAME, _PROGRAM_VERSION))
parser.add_argument('inputfile',
                    help='Input FASTQ file, compression supported')
parser.add_argument('-f', '--filtered',
                    help='Output FASTQ file containing filtered (positive) reads')
parser.add_argument('-u', '--unfiltered',
                    help='Output FASTQ file containing unfiltered (negative) reads')
parser_required_named.add_argument('-i', '--index', required=True,
                                   help='Sequence index to filter for')
parser.add_argument('-m', '--mismatches', default=0, type=int,
                    help='Maximum number of mismatches to tolerate')
parser.add_argument('-c', '--compressed', action='store_true',
                    help='Compress output files (note: file extension not modified)')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Show verbose output')
args = parser.parse_args()


input_path = args.inputfile
out_filtered_path = args.filtered
out_unfiltered_path = args.unfiltered
filter_seq_index = args.index
max_mismatches = args.mismatches
verbose = args.verbose

# HELPER FUNCTIONS
# calculates number of mismatches between s1 and s2
def sum_mismatches(s1, s2):
    if (s1==s2): return 0 # rely on interning to speed up this comparison
    else: return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# writes a fastq entry to handle, avoiding unneeded string formating/concatenation
def write_fastq_entry(handle, title, sequence, quality):
    handle.write('@')
    handle.write(title)
    handle.write('\n')
    handle.write(sequence)
    handle.write('\n+\n')
    handle.write(quality)
    handle.write('\n')


# openers for input and output files
if mimetypes.guess_type(input_path)[1] == 'gzip':
    input_opener = functools.partial(gzip.open, mode='rt')
else:
    input_opener = functools.partial(open, mode='r')

if args.compressed:
    output_opener = functools.partial(gzip.open, mode='wt')
else:
    output_opener = functools.partial(open, mode='w')

# PROCESSING
def main():
    total_reads = 0
    filtered_reads = 0
    unfiltered_reads = 0
    with input_opener(input_path) as input_handle:
        if out_filtered_path:
            filtered_handle = output_opener(out_filtered_path)
        else:
            filtered_handle = None
        if out_unfiltered_path:
            unfiltered_handle = output_opener(out_unfiltered_path)
        else:
            unfiltered_handle = None
        for title, sequence, quality in FastqGeneralIterator(input_handle):
            # Illimina sequence identifier in FASTQ files:
            # see http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
            # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
            # For the Undetermined FASTQ files only, the sequence observed in the index read is written to the FASTQ header in place of the sample number. This information can be useful for troubleshooting demultiplexing.
            # grab the sequence index, use simplistic method for efficiency
            entry_seq_index = title.rsplit(":", 1)[1]
            n_mismatches = sum_mismatches(entry_seq_index, filter_seq_index)
            if verbose:
                print("{} -> index {} -> {} mismatches".format(title,
                                                               entry_seq_index, n_mismatches))
            if n_mismatches <= max_mismatches:
                filtered_reads += 1
                if filtered_handle:
                    write_fastq_entry(filtered_handle, title, sequence, quality)
            else:
                unfiltered_reads += 1
                if unfiltered_handle:
                    write_fastq_entry(unfiltered_handle, title, sequence, quality)
            total_reads += 1

    # OUTPUT
    print("Total reads: {}".format(total_reads))
    print("Filtered reads: {}".format(filtered_reads))
    print("Unfiltered reads: {}".format(unfiltered_reads))

if __name__ == '__main__':
    main()
