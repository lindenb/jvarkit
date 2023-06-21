# SamToJson

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert a SAM input to JSON


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar sam2json  [options] Files

Usage: sam2json [options] Files
  Options:
    -atts, --atts
      do not print attributes
      Default: false
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -cigar, --cigar
      expand cigar
      Default: false
    -flag, --flag
      expand SAm Flags
      Default: false
    -H, --header
      don't print SAM HEADER
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -name, --name
      do not print read name
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * json



## Creation Date

20210402

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamToJson.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamToJson.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam2json** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jvarkit.jar sam2json src/test/resources/toy.bam | python -m json.tool
[
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            },
            {
                "name": "XX",
                "value": [
                    12561,
                    2,
                    20,
                    112
                ]
            }
        ],
        "cigar": "8M4I4M1D3M",
        "flag": 163,
        "len": 39,
        "mapq": 30,
        "matepos": 37,
        "materef": "ref",
        "name": "r001",
        "pos": 7,
        "qualities": "*",
        "ref": "ref",
        "sequence": "TTAGATAAAGAGGATACTG"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "1S2I6M1P1I1P1I4M2I",
        "flag": 0,
        "mapq": 30,
        "name": "r002",
        "pos": 9,
        "qualities": "*",
        "ref": "ref",
        "sequence": "AAAAGATAAGGGATAAA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "5H6M",
        "flag": 0,
        "mapq": 30,
        "name": "r003",
        "pos": 9,
        "qualities": "*",
        "ref": "ref",
        "sequence": "AGCTAA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "6M14N1I5M",
        "flag": 0,
        "mapq": 30,
        "name": "r004",
        "pos": 16,
        "qualities": "*",
        "ref": "ref",
        "sequence": "ATAGCTCTCAGC"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "6H5M",
        "flag": 16,
        "mapq": 30,
        "name": "r003",
        "pos": 29,
        "qualities": "*",
        "ref": "ref",
        "sequence": "TAGGC"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "9M",
        "flag": 83,
        "len": -39,
        "mapq": 30,
        "matepos": 7,
        "materef": "ref",
        "name": "r001",
        "pos": 37,
        "qualities": "*",
        "ref": "ref",
        "sequence": "CAGCGCCAT"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "20M",
        "flag": 0,
        "mapq": 30,
        "name": "x1",
        "pos": 1,
        "qualities": "????????????????????",
        "ref": "ref2",
        "sequence": "AGGTTTTATAAAACAAATAA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "21M",
        "flag": 0,
        "mapq": 30,
        "name": "x2",
        "pos": 2,
        "qualities": "?????????????????????",
        "ref": "ref2",
        "sequence": "GGTTTTATAAAACAAATAATT"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "9M4I13M",
        "flag": 0,
        "mapq": 30,
        "name": "x3",
        "pos": 6,
        "qualities": "??????????????????????????",
        "ref": "ref2",
        "sequence": "TTATAAAACAAATAATTAAGTCTACA"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "25M",
        "flag": 0,
        "mapq": 30,
        "name": "x4",
        "pos": 10,
        "qualities": "?????????????????????????",
        "ref": "ref2",
        "sequence": "CAAATAATTAAGTCTACAGAGCAAC"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "24M",
        "flag": 0,
        "mapq": 30,
        "name": "x5",
        "pos": 12,
        "qualities": "????????????????????????",
        "ref": "ref2",
        "sequence": "AATAATTAAGTCTACAGAGCAACT"
    },
    {
        "atts": [
            {
                "name": "RG",
                "value": "gid1"
            }
        ],
        "cigar": "23M",
        "flag": 0,
        "mapq": 30,
        "name": "x6",
        "pos": 14,
        "qualities": "???????????????????????",
        "ref": "ref2",
        "sequence": "TAATTAAGTCTACAGAGCAACTA"
    }
]
```


