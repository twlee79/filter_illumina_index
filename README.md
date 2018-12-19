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

For information on Illumina sequence identifiers in FASTQ files, see: http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm

### Usage details

```
usage: filter_illumina_index [-h] [--version] [-f FILTERED] [-u UNFILTERED] -i
                             INDEX [-m MISMATCHES] [-c] [-v]
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
                        Maximum number of mismatches to accept (default: 0)
  -c, --compressed      Compress output files (note: file extension not
                        modified) (default: False)
  -v, --verbose         Show verbose output (default: False)

required named arguments:
  -i INDEX, --index INDEX
                        Sequence index to filter for (default: None)
```

### Example usage

The directory `srv` contains example reads in FASTQ and compressed FASTQ format with index `GATCGTGT` and one read with a mismatch.

To test, run:

`filter_illumina_index srv/example_reads.fastq --index GATCGTGT --filtered var/filtered_reads.fastq --unfiltered var/unfiltered_reads.fastq`

This will process `srv/example_reads.fastq`, matching to index `GATCGTGT` with no mismatches allowed (default). Reads matching this index will be saved to `var/filtered_reads.fastq` and those not matching this index will be saved to  `var/unfiltered_reads.fastq`. In addition, the following output will be displayed:

```
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
```

---

### Additional details

* Author:       Tet Woo Lee
* Copyright:    © 2018 Tet Woo Lee
* Licence:      GPLv3
* Dependencies: Biopython, tested on v1.72

### Change log

version 1.0.2 2018-12-19
: Shows statistics on number of mismatches found

version 1.0.1 2018-12-19
: Speed up number of mismatches calculation

version 1.0 2018-12-14
: Minor updates for PyPi and conda packaging

version 1.0.dev1 2018-12-13
: First working version
