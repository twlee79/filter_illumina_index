# filter_illumina_index
## Filter a Illumina FASTQ file based on index sequence.

Reads a Illumina FASTQ file and compares the sequence index in the
`sample number` position of the sequence identifier to a supplied sequence
index. Entries that match the sequence index are filtered into the *filtered
file* (if any) and entries that don't match are filtered into the *unfiltered
file* (if any). Displays the count of total, filtered and unfiltered reads,
as well as the number of mismatches found across all reads. Matching tolerating
a certain number of mismatches (`-m` parameter), and gzip compression for input
(detected on the basis of file extension) and output (specified using `-c`
parameter) are supported.

Specifying an empty index, (`-i ""` or `-i` with no argumenent) enables
'passthrough' mode where all reads are directed to the output filtered file with
no processing. Passthrough mode is useful if this program is part of a workflow
that needs to be adapted to files that do not have a valid Illumina index, as it
allows all processing of this program to be skipped.

For information on Illumina sequence identifiers in FASTQ files, see:
http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm

### Usage details

```
usage: filter_illumina_index [-h] [--version] [-f FILTERED] [-u UNFILTERED] -i
                             [INDEX] [-m MISMATCHES] [-c] [-v]
                             inputfile

positional arguments:
  inputfile             Input FASTQ file, compression supported

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -f FILTERED, --filtered FILTERED
                        Output FASTQ file containing filtered (positive) reads
                        (default: None)
  -u UNFILTERED, --unfiltered UNFILTERED
                        Output FASTQ file containing unfiltered (negative)
                        reads (default: None)
  -m MISMATCHES, --mismatches MISMATCHES
                        Maximum number of mismatches to tolerate (default: 0)
  -c, --compressed      Compress output files (note: file extension not
                        modified) (default: False)
  -v, --verbose         Show verbose output (default: False)

required named arguments:
  -i [INDEX], --index [INDEX]
                        Sequence index to filter for; if empty (i.e. no
                        argument or "") then program will run in "passthrough"
                        mode with all reads directed to filtered file with no
                        processing (default: None)
```

### Example usage

The directory `` contains example reads in FASTQ and compressed FASTQ format with index `GATCGTGT` and one read with a mismatch.

To test, run:

`filter_illumina_index usr/example_reads.fastq --index GATCGTGT --filtered var/filtered_reads.fastq --unfiltered var/unfiltered_reads.fastq`

python filter_illumina_index/filter_illumina_index.py usr/example_reads.fastq --index GATCGTGT 
filter_illumina_index 1.0.4.dev1
Input file: usr/example_reads.fastq
Filtering for sequence index: GATCGTGT
Max mismatches tolerated: 0
Output filtered file: None
Output unfiltered file: None
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

This will process `srv/example_reads.fastq`, matching to index `GATCGTGT` with
no mismatches allowed (default). Reads matching this index will be saved to
`var/filtered_reads.fastq` and those not matching this index will be saved to
`var/unfiltered_reads.fastq`. In addition, the following output will be
displayed:

```
filter_illumina_index 1.0.4.dev1
Input file: usr/example_reads.fastq
Filtering for sequence index: GATCGTGT
Max mismatches tolerated: 0
Output filtered file: None
Output unfiltered file: None
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

The number of mismatches is simply the number of characters that differ between
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
  xopen, tested with v0.8.4


### Change log

version 1.0.3.post2 2020-01-04  
Improved output
  - Add version number to output
  - Show parameters in output
  - Allow no argument for passthrough mode

version 1.0.3.post1 2020-01-04  
Minor bugfix
  - Bugfix: Bump version number in script

version 1.0.3 2020-01-04  
Added `passthrough` mode with empty index

version 1.0.2 2018-12-19  
Shows statistics on number of mismatches found

version 1.0.1 2018-12-19  
Speed up number of mismatches calculation

version 1.0 2018-12-14  
Minor updates for PyPi and conda packaging

version 1.0.dev1 2018-12-13  
First working version
