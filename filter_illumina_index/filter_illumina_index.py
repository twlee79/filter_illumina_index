#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import functools
import mimetypes
import gzip
import sys
from functools import partial

import dnaio
import xopen

"""
Filter a Illumina FASTQ file based on index sequence.

Reads a Illumina FASTQ file and compares the sequence index in the
`sample number` position of the sequence identifier to a supplied sequence
index. Entries that match the sequence index are filtered into the *filtered
file* (if any) and entries that don't match are filtered into the *unfiltered
file* (if any). Displays the count of total, filtered and unfiltered reads,
as well as the number of mismatches found across all reads. Matching tolerating
a certain number of mismatches (`-m` parameter), and gzip compression for input
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
# Dependencies: dnaio, tested with v0.4.1
#               xopen, tested with v0.9.0
# -------------------------------------------------------------------------------

_PROGRAM_VERSION = '1.0.4.dev2'
# -------------------------------------------------------------------------------
# ### Change log
#
# version 1.0.4 2020-01-04
# Speed up and algorithm changes
#   - Switch to `dnaio` over Biopython to improve speed (>3x faster + multi-
#     threading support for compression)
#   - Change mismatch calculation algorithm, now includes any characters
#     missing in filter-by index or read index
#
# version 1.0.3.post2 2020-01-04
# Improved output
#   - Add version number to output
#   - Show parameters in output
#   - Allow no argument for passthrough mode
#
# version 1.0.3.post1 2020-01-04
# Minor bugfix
#   - Bugfix: Bump version number in script
#
# version 1.0.3 2020-01-04
# Added `passthrough` mode with empty index
#
# version 1.0.2 2018-12-19
# Shows statistics on number of mismatches found
#
# version 1.0.1 2018-12-19
# Speed up number of mismatches calculation
#
# version 1.0 2018-12-14
# Minor updates for PyPi and conda packaging
#
# version 1.0.dev1 2018-12-13
# First working version
#
# -------------------------------------------------------------------------------

# TODO: tests

# INITIALISATION

def main(argv = None, return_result = False):
        # return_result = True will return summary of output to caller (for testing)
    if argv is None: argv = sys.argv[1:] # if parameters not provided, use sys.argv

    _PROGRAM_NAME_VERSION = '{} {}'.format(_PROGRAM_NAME, _PROGRAM_VERSION)
    parser = argparse.ArgumentParser(prog=_PROGRAM_NAME,
                                     description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_required_named = parser.add_argument_group('required named arguments')
    parser.add_argument('--version', action='version',
                        version=_PROGRAM_NAME_VERSION)
    parser.add_argument('inputfile',
                        help='Input FASTQ file, compression (`.gz`, `.bz2` and '
                        '`.xz`) supported')
    parser.add_argument('-f', '--filtered',
                        help='Output FASTQ file containing filtered (positive) reads; '
                        'compression detected by extension')
    parser.add_argument('-u', '--unfiltered',
                        help='Output FASTQ file containing unfiltered (negative) reads; '
                        'compression detected by extension')
    parser_required_named.add_argument('-i', '--index', required=True, const='', nargs='?',
                                       help='Sequence index to filter for; if empty '
                                       '(i.e. no argument or "") then program will '
                                       'run in "passthrough" mode with all reads '
                                       'directed to filtered file with no processing')
    parser.add_argument('-m', '--mismatches', default=0, type=int,
                        help='Maximum number of mismatches to tolerate')
    parser.add_argument('-t', '--threads', default=1, type=int,
                        help='Number of threads to pass to `xopen` for each '
                        'open file; use 0 to turn off `pigz` use and rely '
                        'on `gzip.open` so no extra threads spawned.')
    parser.add_argument('-l', '--compresslevel', default=6, type=int, choices=range(1,10),
                        help='Compression level for writing gzip files; '
                        'ignored if gzip compression not used')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='Increase logging verbosity, available levels 1 to 2 '
                             'with `-v` to `-vv`')
    args = parser.parse_args(argv)


    input_path = args.inputfile
    out_filtered_path = args.filtered
    out_unfiltered_path = args.unfiltered
    filter_seq_index = args.index
    max_tolerated_mismatches = args.mismatches
    threads = args.threads
    compresslevel = args.compresslevel
    verbose = args.verbose

    if filter_seq_index == '':
        passthrough_mode = True
        max_tolerated_mismatches = float('NaN')
    else:
        passthrough_mode = False

    # HELPER FUNCTIONS
    print(_PROGRAM_NAME_VERSION)
    print("Input file: {}".format(input_path))
    print("Filtering for sequence index: {}{}".format(filter_seq_index,
        "(passthrough mode)" if passthrough_mode else ""))
    print("Max mismatches tolerated: {}".format(max_tolerated_mismatches))
    print("Output filtered file: {}".format(out_filtered_path))
    print("Output unfiltered file: {}".format(out_unfiltered_path))
    if verbose>=1:
        print("Using {} threads per open file".format(threads))
        print("Compression level: {}".format(compresslevel))
        print("Showing verbose level {} logging".format(verbose))


    xopen_xthreads = functools.partial(xopen.xopen, threads=threads,
                                       compresslevel = compresslevel)

    # PROCESSING
    total_reads = 0
    filtered_reads = 0
    unfiltered_reads = 0
    max_tracked_mismatches = len(filter_seq_index)
    cumul_n_mismatches = [0 for i in range(max_tracked_mismatches + 2)]
    # array for tracking number of mismatches (0-indexlen,>indexlen+1)
    # >indexlen+1 is required because we define mismatch to include extra
    # characters from desired index OR read index, which may be longer
    if passthrough_mode: filtered = True
    with dnaio.open(input_path, mode='r', opener=xopen_xthreads) as input_fastq:
        if out_filtered_path:
            filtered_fastq = dnaio.open(out_filtered_path, mode='w', opener=xopen_xthreads)
        else:
            filtered_fastq = None
        if out_unfiltered_path:
            unfiltered_fastq = dnaio.open(out_unfiltered_path, mode='w', opener=xopen_xthreads)
        else:
            unfiltered_fastq = None
        for record in input_fastq:
            total_reads += 1
            if passthrough_mode:
                # always filtered
                if verbose>=2:
                    seqid = record.name
                    print("{} (passthrough-mode) (filtered)".format(seqid))

            else:
                # Illimina sequence identifier in FASTQ files:
                # see http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
                # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
                # For the Undetermined FASTQ files only, the sequence observed in the index read is written to the FASTQ header in place of the sample number. This information can be useful for troubleshooting demultiplexing.
                # grab the sequence index, use simplistic method of value after last : for efficiency
                seqid = record.name
                entry_seq_index = seqid[seqid.rfind(':')+1:]
                n_mismatches = abs(len(entry_seq_index)-len(filter_seq_index))
                if entry_seq_index!=filter_seq_index:
                    for c1,c2 in zip(entry_seq_index,filter_seq_index):
                        if c1!=c2: n_mismatches+=1
                filtered = (n_mismatches <= max_tolerated_mismatches)
                if verbose>=2:
                    print("{} -> index {} -> {} mismatches ({})".format(seqid,
                        entry_seq_index, n_mismatches,
                        'filtered' if filtered else 'unfiltered'))
                if n_mismatches>=max_tracked_mismatches:
                    n_mismatches = max_tracked_mismatches+1
                cumul_n_mismatches[n_mismatches] += 1
            if filtered:
                filtered_reads += 1
                if filtered_fastq: filtered_fastq.write(record)
            else:
                unfiltered_reads += 1
                if unfiltered_fastq: unfiltered_fastq.write(record)
        if filtered_fastq: filtered_fastq.close()
        if unfiltered_fastq: unfiltered_fastq.close()


    # OUTPUT
    print("Total reads: {}".format(total_reads))
    if passthrough_mode:
        print("Filtered reads: {} (passthrough-mode)".format(filtered_reads))
    else:
        print("Filtered reads: {}".format(filtered_reads))
        print("Unfiltered reads: {}".format(unfiltered_reads))
        for n_mismatches, cumul_mismatches in enumerate(cumul_n_mismatches):
            print(" Reads with {}{} mismatches: {}".format(
                '' if n_mismatches<= max_tracked_mismatches else '>=',
                n_mismatches, cumul_mismatches))

    if return_result: 
        results = {
            "total" : total_reads,
            "filtered": filtered_reads,
            "unfiltered": unfiltered_reads,
            "mismatches": {n_mismatches : cumul_mismatches 
                          for n_mismatches,cumul_mismatches in enumerate(cumul_n_mismatches)}
        }
        return(results)

if __name__ == '__main__':
    main()
