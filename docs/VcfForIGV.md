# VcfForIGV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Prepare IGV sessions file and json for nextflow/igv_report


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcf4igv  [options] Files

Usage: vcf4igv [options] Files
  Options:
    --count
      For Case-controls (at least one 'case' and one 'control' in the 
      samplesheet' ) 8 comma integers representing the maximum number of bam 
      to display: case-HOM_REF,case-HET,case-HOM_VAR,case-NO_CALL,ctrl-HOM_REF,ctrl-HET,ctrl-HOM_VAR,ctr-NO_CALL 
      . If there is no case  4 comma separated 
      integers:HOM_REF,HET,HOM_VAR,NO_CALL (other fields are ignored)
      Default: 1,5,1,0,1,5,1,0
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
  * --samplesheet
      samplesheet. TSV or CSV file with the following header: 
      'bam,sample,status'. 'bam' is required, 'status' must be 'case' or 
      'control' 
    --session-dir
      if defined, save the IGV session XML files in that directory
    --version
      print version and exit

```


## Keywords

 * vcf
 * igv
 * json



## Creation Date

20260506

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf4igv/VcfForIGV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf4igv/VcfForIGV.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcf4igv/VcfForIGVTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcf4igv/VcfForIGVTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf4igv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example
 
 ```
 java -jar dist/jvarkit.jar vcf4igv -R src/test/resources/rotavirus_rf.fa --samplesheet jeter.csv  src/test/resources/rotavirus_rf.vcf.gz | python3 -m json.tool

[
    {
        "fasta": "/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa",
        "interval": "RF01:970-970",
        "contig": "RF01",
        "chromosome": "RF01",
        "start": 970,
        "end": 970,
        "length": 1,
        "ref": "A",
        "alt": "C",
        "description": "RF01:970 CASE/HET:N=0.  CASE/HOM_VAR:N=1 displayed here: S5.  CASE/HOM_REF:N=2 displayed here: S1.  CTRL/HET:N=0.  CTRL/HOM_VAR=0.  CTRL/HOM_REF:N=2 displayed here: S2. ",
        "bams": [
            {
                "sample": "S5",
                "bam": "src/test/resources/S5.bam",
                "status": "case"
            },
            {
                "sample": "S1",
                "bam": "src/test/resources/S1.bam",
                "status": "case"
            },
            {
                "sample": "S2",
                "bam": "src/test/resources/S2.bam",
                "status": "control"
            }
        ]
    },
(...)
```


