#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
filter_illumina_index Unit Tests
========================

Performs the following tests:
    - Loading the script with no options
    - Loading the script with --version
    - Running the script with various options using a test input file. The
      output of the script is compared line-by-line to reference output files.

To run tests, use one of these commands:
|   python -m unittest discover
|   python -m filter_illumina_index.tests.test_filter_illumina_index

The reference output files can be (re)generated with the following command.
Suggest manually checking these output files (using git to easily see changes)
prior to packaging these in a new version:
|   python -m filter_illumina_index.tests.test_filter_illumina_index --generate

"""

# Currently does the following tests:
#   Script loading and argparse exitcodes
#   Reading fastq and fastq.gz
#   Mismatch counting by comparing summaries (compared to expected)
#     With different index
#   Writing fastq and fastq.gz to filtered/unfiltered (compared to expected)
#     With different index and diff # tolerated mismatches

import unittest
import sys
import functools
import json
import itertools
import gzip

from filter_illumina_index.filter_illumina_index import main as filter_illumina_index_main2
filter_illumina_index_main = functools.partial(filter_illumina_index_main2, return_result = True)


tests_root = 'filter_illumina_index/tests/data/'
tests_results_root = 'filter_illumina_index/tests/data/results/'
tests_output_root = 'filter_illumina_index/tests/var/'
input_test_file_fastq = tests_root + 'test_reads_GATCGTGT.fastq'
input_test_file_fastq_gz = tests_root + 'test_reads_GATCGTGT.fastq.gz'
input_test_file_diffbarcodes = tests_root + 'test_reads_diffbarcodes.fastq'


test_set_exitcodes = [
    # tuples of ([options], expected return code)
    # only for sys.exit by argparse
    ([], 2),
    (['--version'], 0),
    ([input_test_file_fastq], 2),
    (['--index','GATCGTGT'], 2),
]

test_sets_vs_summary = [
    # tuples of ([options], expected output file,  [generate])
    # generate is optional, if False, then the results will not be generated
    # for this test set

    # READ TESTS
    # test read fastq to expected summary, with dnaio -t 1 (default) or -t 0
    ([input_test_file_fastq, '--index','GATCGTGT'],'test_reads_GATCGTGT_results.json'),
    ([input_test_file_fastq, '--index','GATCGTGT','-t','0'],'test_reads_GATCGTGT_results.json',False),
    # test read fastq.gz to expected summary, with dnaio -t 1 (default) or -t 0
    ([input_test_file_fastq_gz, '--index','GATCGTGT'],'test_reads_GATCGTGT_results.json',False),
        # this gives a ResourceWarning for unclosed files in xopen v0.8.4, fixed in v0.9.0
    ([input_test_file_fastq_gz, '--index','GATCGTGT','-t','0'],'test_reads_GATCGTGT_results.json',False),

    # MISMATCH NUMBER TESTS
    # test of mismatch calculation, see examples in README
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv'],'test_reads_diffbarcodes_results_TGACCAAT.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv'],'test_reads_diffbarcodes_results_NNNCCAAT.json'),

    # MISMATCH FILTER/UNFILTER TESTS
    # test number of filtered/unfiltered with diff mismatch tolerances
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 1'],'test_reads_diffbarcodes_results_TGACCAAT_m1.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 1'],'test_reads_diffbarcodes_results_NNNCCAAT_m1.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 2'],'test_reads_diffbarcodes_results_TGACCAAT_m2.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 2'],'test_reads_diffbarcodes_results_NNNCCAAT_m2.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 3'],'test_reads_diffbarcodes_results_TGACCAAT_m3.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 3'],'test_reads_diffbarcodes_results_NNNCCAAT_m3.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 4'],'test_reads_diffbarcodes_results_TGACCAAT_m4.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 4'],'test_reads_diffbarcodes_results_NNNCCAAT_m4.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 5'],'test_reads_diffbarcodes_results_TGACCAAT_m5.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 5'],'test_reads_diffbarcodes_results_NNNCCAAT_m5.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 6'],'test_reads_diffbarcodes_results_TGACCAAT_m6.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 6'],'test_reads_diffbarcodes_results_NNNCCAAT_m6.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 7'],'test_reads_diffbarcodes_results_TGACCAAT_m7.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 7'],'test_reads_diffbarcodes_results_NNNCCAAT_m7.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 8'],'test_reads_diffbarcodes_results_TGACCAAT_m8.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 8'],'test_reads_diffbarcodes_results_NNNCCAAT_m8.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 10'],'test_reads_diffbarcodes_results_TGACCAAT_m10.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 10'],'test_reads_diffbarcodes_results_NNNCCAAT_m10.json'),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-vv','-m 20'],'test_reads_diffbarcodes_results_TGACCAAT_m20.json'),
    ([input_test_file_diffbarcodes, '--index','NNNCCAAT','-vv','-m 20'],'test_reads_diffbarcodes_results_NNNCCAAT_m20.json'),
]

test_sets_vs_output = [
    # tuples of ([options], (expected filtered file, expected unfiltered file), [generate])
    # either filtered/unfiltered can be None in which not checked against
    # end with .gz to automatically test compressed outputs 
    # generate is optional, if False, then the results will not be generated
    # for this test set
    ([input_test_file_fastq, '--index','GATCGTGT',],
     ('test_reads_GATCGTGT_filtered.fastq','test_reads_GATCGTGT_unfiltered.fastq')
    ),
    ([input_test_file_fastq, '--index','GATCGTGT',],
     (None,'test_reads_GATCGTGT_unfiltered.fastq'), False
    ),
    ([input_test_file_fastq, '--index','GATCGTGT',],
     ('test_reads_GATCGTGT_filtered.fastq', None), False
    ),
    ([input_test_file_fastq, '--index','GATCGTGT',],
     ('test_reads_GATCGTGT_filtered.fastq.gz','test_reads_GATCGTGT_unfiltered.fastq.gz')
    ),
    ([input_test_file_fastq, '--index','GATCGTCT',], # inverted, GATCGTCT is mismatch
     ('test_reads_GATCGTGT_unfiltered.fastq.gz','test_reads_GATCGTGT_filtered.fastq.gz'), False
    ),

    # tests with mismatches allowed
    ([input_test_file_fastq, '--index','GATCGTGT','-m','1'], # allow 1 mismatch
     ('test_reads_GATCGTGT_filtered_m1.fastq.gz','test_reads_GATCGTGT_unfiltered_m1.fastq')
    ),
    ([input_test_file_fastq, '--index','GATCGTCT','-m','1'], # inverted, GATCGTCT is mismatch, but allow 1 mismatch
     ('test_reads_GATCGTGT_filtered_m1.fastq.gz','test_reads_GATCGTGT_unfiltered_m1.fastq'), False
    ),

    # tests with varyinf # mismatches allowed
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','0','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m0.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m0.fastq.gz'), 
    ),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','1','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m1.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m1.fastq.gz'), 
    ),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','2','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m2.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m2.fastq.gz'), 
    ),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','3','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m3.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m3.fastq.gz'), 
    ),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','5','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m5.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m5.fastq.gz'), 
    ),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','8','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m8.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m8.fastq.gz'), 
    ),
    ([input_test_file_diffbarcodes, '--index','TGACCAAT','-m','10','-vv'],
     ('test_reads_diffbarcodes_filtered_TGACCAAT_m10.fastq.gz',
      'test_reads_diffbarcodes_unfiltered_TGACCAAT_m10.fastq.gz'), 
    )

]

class TestFilerIlluminaIndex(unittest.TestCase):
    def test_exitcodes(self):
        for test_set in test_set_exitcodes:
            test_options, test_expected_exitcode = test_set
            flat_test_options = " ".join(test_options)
            with self.subTest(options = flat_test_options,
                              exitcode = test_expected_exitcode):
                print("Testing options {}, expecting exitcode {}:".format(
                    flat_test_options,
                    test_expected_exitcode))
                with self.assertRaises(SystemExit) as cm:
                    filter_illumina_index_main(test_options)
                exitcode = cm.exception.code
                self.assertEqual(exitcode, test_expected_exitcode, 
                    'Unexpected return code: got {}, expected {}'.format(exitcode, test_expected_exitcode))

    def test_results_summary(self):
        # generic tests compared to saved summary
        for test_set in test_sets_vs_summary:
            test_options, test_expected_output_file = test_set[:2]
            flat_test_options = " ".join(test_options)
            test_expected_output_path = tests_results_root + test_expected_output_file
            with self.subTest(options = flat_test_options,
                              reference_file = test_expected_output_file):
                print("Testing options {}, comparing summary to {}:".format(
                    flat_test_options,
                    test_expected_output_file))
                test_result = filter_illumina_index_main(test_options)
                test_json_safe = json.loads(json.dumps(test_result))
                    # convert to json and back to ensure it is comparable to jsonified expected results
                with open(test_expected_output_path,'r') as expected_obj:
                    expected_result = json.load(expected_obj)
                self.assertEqual(test_json_safe, expected_result,'Mismatch results {} vs {}'.format(
                    test_json_safe, expected_result))

    def helper_compare_files(self, path1, path2):
        opener1 = functools.partial(gzip.open, mode='rt') if path1.endswith('gz') else functools.partial(open, mode='r')
        opener2 = functools.partial(gzip.open, mode='rt') if path2.endswith('gz') else functools.partial(open, mode='r')
        with opener1(path1) as handle1:
            with opener2(path2) as handle2:
                for line_index,(line1,line2) in enumerate(
                    itertools.zip_longest(handle1, handle2),1):
                    self.assertEqual(line1.strip('\r\n'),
                        line2.strip('\r\n'),
                        '{} vs {}, line {} mismatch'.format(path1,
                            path2, line_index))

    def test_results_output(self):
        for test_set in test_sets_vs_output:
            test_options, test_expected_output_pair = test_set[:2]
            flat_test_options = " ".join(test_options)
            test_expected_filtered_file, test_expected_unfiltered_file = test_expected_output_pair
            if test_expected_filtered_file is not None:
                test_expected_filtered_path = tests_results_root + test_expected_filtered_file
                test_output_filtered_path = tests_output_root + test_expected_filtered_file
                test_options.extend(['-f',test_output_filtered_path])
            if test_expected_unfiltered_file is not None:
                test_expected_unfiltered_path = tests_results_root + test_expected_unfiltered_file
                test_output_unfiltered_path = tests_output_root + test_expected_unfiltered_file
                test_options.extend(['-u',test_output_unfiltered_path])
            with self.subTest(options = flat_test_options,
                              filtered_file = test_expected_filtered_file,
                              unfiltered_file = test_expected_unfiltered_file,
                              ):
                print("Testing options {}, comparing output to filtered {} / unfiltered{}:".format(
                    flat_test_options,
                    test_expected_filtered_file,
                    test_expected_unfiltered_file))
                filter_illumina_index_main(test_options)
                if test_expected_filtered_file:
                    self.helper_compare_files(test_expected_filtered_path, test_output_filtered_path)
                if test_expected_unfiltered_file:
                    self.helper_compare_files(test_expected_unfiltered_path, test_output_unfiltered_path)


# this function will generate results files, instead of checking against them
# will be run if called with option --generate
def generate_results():
    for test_set in test_sets_vs_summary:
        if len(test_set)>2: # third parameter == False, don't generate
            if test_set[2]  == False:
                continue
        test_options, test_expected_output_file = test_set[:2]
        test_expected_output_path = tests_results_root + test_expected_output_file
        print("Running with options {}, saving to {}:".format(' '.join(test_options), 
            test_expected_output_path))
        results = filter_illumina_index_main(test_options)
        with  open(test_expected_output_path,'w') as test_expected_output_handle:
            json.dump(results, test_expected_output_handle, indent=1)
            
    for test_set in test_sets_vs_output:
        if len(test_set)>2: # third parameter == False, don't generate
            if test_set[2]  == False:
                continue
        test_options, test_expected_output_pair = test_set[:2]
        test_expected_filtered_file, test_expected_unfiltered_file = test_expected_output_pair
        if test_expected_filtered_file is not None:
            test_expected_filtered_path = tests_results_root + test_expected_filtered_file
            test_options.extend(['-f',test_expected_filtered_path])
        if test_expected_unfiltered_file is not None:
            test_expected_unfiltered_path = tests_results_root + test_expected_unfiltered_file
            test_options.extend(['-u',test_expected_unfiltered_path])
        print("Running with options {}:".format(' '.join(test_options)))
        filter_illumina_index_main(test_options)


if __name__ == '__main__':
    if '--generate' in sys.argv:
        do_generate = True
    else:
        do_generate = False
    if do_generate: generate_results()
    else: unittest.main()
