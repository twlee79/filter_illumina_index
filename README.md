# filter_illumina_index
## Filter a Illumina FASTQ file based on index sequence.

Reads a Illumina FASTQ file and compares the sequence index in the
`sample number` position of the sequence identifier to a supplied sequence
index. Entries that match the sequence index are filtered into the *filtered
file* (if any) and entries that don't match are filtered into the *unfiltered
file* (if any). Displays the count of total, filtered and unfiltered reads,
as well as the number of mismatches found across all reads. Matching tolerating
a certain number of mismatches (`-m` parameter), and compression for input
and output are supported (detected on the basis of file extension).

Specifying an empty index, (`-i ""` or `-i` with no argument) enables
'passthrough' mode where all reads are directed to the output filtered file with
no processing. Passthrough mode is useful if this program is part of a workflow
that needs to be adapted to files that do not have a valid Illumina index, as it
allows all processing of this program to be skipped.

### Illumina index

For information on Illumina sequence identifiers in FASTQ files, see [FASTQ 
Files from Illumina](https://help.basespace.illumina.com/articles/descriptive/fastq-files/).
This includes the following excerpt:

```
For the Undetermined FASTQ files only, the sequence observed in the index read 
is written to the FASTQ header in place of the sample number. This information 
can be useful for troubleshooting demultiplexing.

@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
```

This script assumes the `<sample number>` value is the Illumina index. To
ensure rapid processing, the value following the final colon(`:`) is used
without confirming the sequence identifier conforms to the above format, or
whether the index is actually a nucleotide sequence.

### Usage details

```
usage: filter_illumina_index [-h] [--version] [-f FILTERED] [-u UNFILTERED] -i
                             [INDEX] [-m MISMATCHES] [-t THREADS]
                             [-l {1,2,3,4,5,6,7,8,9}] [-v]
                             inputfile

positional arguments:
  inputfile             Input FASTQ file, compression (`.gz`, `.bz2` and
                        `.xz`) supported

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -f FILTERED, --filtered FILTERED
                        Output FASTQ file containing filtered (positive)
                        reads; compression detected by extension (default:
                        None)
  -u UNFILTERED, --unfiltered UNFILTERED
                        Output FASTQ file containing unfiltered (negative)
                        reads; compression detected by extension (default:
                        None)
  -m MISMATCHES, --mismatches MISMATCHES
                        Maximum number of mismatches to tolerate (default: 0)
  -t THREADS, --threads THREADS
                        Number of threads to pass to `xopen` for each open
                        file; use 0 to turn off `pigz` use and rely on
                        `gzip.open` so no extra threads spawned. (default: 1)
  -l {1,2,3,4,5,6,7,8,9}, --compresslevel {1,2,3,4,5,6,7,8,9}
                        Compression level for writing gzip files; ignored if
                        gzip compression not used (default: 6)
  -v, --verbose         Increase logging verbosity, available levels 1 to 2
                        with `-v` to `-vv` (default: 0)

required named arguments:
  -i [INDEX], --index [INDEX]
                        Sequence index to filter for; if empty (i.e. no
                        argument or "") then program will run in "passthrough"
                        mode with all reads directed to filtered file with no
                        processing (default: None)
```

### Example usage

The directory `filter_illumina_index/tests/data/` (within the location where
the package is installed) contains various test files, for example 
`test_reads_GATCGTGT.fastq` with 29 reads containing barcode `GATCGTGT` and
one read with a 1-nt mismatch of `AATCGTGT`. Running this program with the 
following command:

`filter_illumina_index filter_illumina_index/tests/data/test_reads_GATCGTGT.fastq --index GATCGTGT --filtered /tmp/filtered_reads.fastq --unfiltered /tmp/unfiltered_reads.fastq`

will process that input file with no mismatches allowed (default) for the index
`GATCGTGT`. The output files are `/tmp/filtered_reads.fastq` with 29 reads 
containing the `GATCGTGT` barcode, and `/tmp/unfiltered_reads.fastq` with 1 read
containing the mismatched-barcode `AATCGTGT`. The following log to stdout will 
be displayed:

```
filter_illumina_index 1.0.4
Input file: filter_illumina_index/tests/data/test_reads_GATCGTGT.fastq
Filtering for sequence index: GATCGTGT
Max mismatches tolerated: 0
Output filtered file: /tmp/filtered_reads.fastq
Output unfiltered file: /tmp/unfiltered_reads.fastq
Total reads: 30
Filtered reads: 29
Unfiltered reads: 1
 Reads with 0 mismatches: 29
 Reads with 1 mismatches: 1
 Reads with 2 mismatches: 0
 Reads with 3 mismatches: 0
 Reads with 4 mismatches: 0
 Reads with 5 mismatches: 0
 Reads with 6 mismatches: 0
 Reads with 7 mismatches: 0
 Reads with 8 mismatches: 0
 Reads with >=9 mismatches: 0
```

### Algorithm details

The barcode is read from the sequence number position of the sequence identifier
using a very simplistic method to optimise performance: the field in the
sequence identifier following the last colon (`:`) character is used, e.g.

```
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1 1:N:0:TGACCAAT
                                          ^         : last colon
                                           ^^^^^^^^ : taken
                                           TGACCAAT : barcode
```

If there is no colon in the sequence identifier, an exception is produced as
no barcode can be found, except in passthrough mode (which does no actual
processing and redirects all reads to the filtered file). The number of 
mismatches is simply the number of characters that differ between
the provided index sequence and the barcode from the read. Any missing
characters in either the provided index or the barcode are counted as
mismatches. The comparison is performed at the start of the each set of
characters with no provision for wildcards, insertions or deletions, e.g.

```
TGACCAAT index
TGACCAAT read barcode
         0 mismatches

TGACCAAT index
AGACCAAT read barcode
^        1 mismatch

TGACCAAT index
GACCAAT  read barcode
^^^ ^ ^^ 6 mismatches

TGACCAAT index
TGACCAAAA read barcode
       ^^ 2 mismatches

TGACCAAT index
NNNCCAAT read barcode
^^^      3 mismatches

NNNCCAAT index
TGACCAAT read barcode
^^^      3 mismatches

NNNCCAAT index
NNNCCAAT read barcode
         0 mismatches

TGACCAAT index
TGACCAATTGACCAATTGACCAAT read barcode
        ^^^^^^^^^^^^^^^^ 16 mismatches
```

Where there is a greater number of characters in the read barcode than the
provided index, the number of mismatches is summarised as `>=index length+1`,
i.e. the final entry above will be counted as >=9 mismatches.

### File reading/writing and threading

This script uses the `dnaio` and `xopen` packages for reading/writing FASTQ
files with compression support. The `xopen` package used for reading/writing
compressed files spawns `pigz` processes to speed-up processing. The `--threads`
parameter indicates the number of threads passed to the `xopen()` function. With
`--threads 1` there is still a small amount of multithreading as a different
`pigz` process is spawned for each open file and up to 3 files are open at once
(one input and two output), although this is largely throttled by the sequential
nature of the script processing. To prevent these `pigz` processes being
spawned, `--threads 0` can be used, which causes a fallback to `gzip.open` at
an additional performance cost (it is slower than `pigz`).

`xopen` supports automatic compression according to the file extension, with
`.gz`, `.bz2` and `.xz` supported. Only the first of these has been tested
with this script but the others are expected to work without issue.

---

### Additional details

* Author:       Tet Woo Lee
* Copyright:    © 2018-2020 Tet Woo Lee
* Licence:      GPLv3
* Dependencies:  
  dnaio, tested with v0.4.1  
  xopen, tested with v0.9.0


### Change log

version 1.0.4 2020-04-11  
Speed up and algorithm changes
  - Switch to `dnaio` over Biopython to improve speed (>3x faster + multi-
    threading support for compression)
  - Change mismatch calculation algorithm, now includes any characters
    missing in filter-by index or read index
  - Exception if no barcode detected outside of passthrough mode
  - Add unit tests

version 1.0.3.post2 2020-04-01  
Improved output
  - Add version number to output
  - Show parameters in output
  - Allow no argument for passthrough mode

version 1.0.3.post1 2020-04-01  
Minor bugfix
  - Bugfix: Bump version number in script

version 1.0.3 2020-04-01  
Added `passthrough` mode with empty index

version 1.0.2 2018-12-19  
Shows statistics on number of mismatches found

version 1.0.1 2018-12-19  
Speed up number of mismatches calculation

version 1.0 2018-12-14  
Minor updates for PyPi and conda packaging

version 1.0.dev1 2018-12-13  
First working version

