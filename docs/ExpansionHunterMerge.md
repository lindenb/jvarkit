# ExpansionHunterMerge

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge Vcf from ExpansionHunter.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar expansionhuntermerge  [options] Files

Usage: expansionhuntermerge [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --cases
      File or comma-separated list of control samples
    --controls
      File or comma-separated list of control samples
    --factor
      multiple median/mean value of controls by 'factor'. if median value=100, 
      then we count case having a size greater than 100*factor for the burden 
      test 
      Default: 1.0
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    -p, --percentile
      percentile to use 'average' or 'median'
      Default: median
    --skip-filtered
      Skip filtered variants
      Default: false
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --types
      types of call to consider using FORMAT/SO. There can be a bias of size 
      depending of the nature of the call: SPANNING|FLANKING|INREPEAT
      Default: SPANNING,FLANKING,INREPEAT
    --version
      print version and exit

```


## Keywords

 * vcf
 * merge
 * ExpansionHunter



## Creation Date

20210210

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/expansionhunter/ExpansionHunterMerge.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/expansionhunter/ExpansionHunterMerge.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **expansionhuntermerge** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
# Input
 
Input is a list of indexed vcf files or one file with the '.list' suffix containing the path to the vcfs
 


