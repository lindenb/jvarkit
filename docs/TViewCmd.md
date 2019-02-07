# TViewCmd

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

equivalent of samtools tview


## Usage

```
Usage: tview [options] Files
  Options:
    --clip
      Show clip
      Default: false
    --coverage, --depth
      Number of rows for coverage (hide:<=0)
      Default: 10
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hideBases
      Hide bases
      Default: false
    --hideHomRef
      Hide HOM_REF variations
      Default: false
    --hideNoCall
      Hide NO_CALL variations
      Default: false
    --insert
      Show insertions
      Default: false
    -layout, --layout
      Layout reads
      Default: pileup
      Possible Values: [pileup, name]
    -left, --leftmargin
      left margin width
      Default: 15
    -maxrows, --maxrowss
      maximum number of rows per read group. -1 == all
      Default: -1
    --noconsensus
      Hide Consensus line
      Default: false
    --nodefaultinterval
      if no interval was provided, don't try to create a default one.
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --format, --outputformat
      Output format
      Default: tty
      Possible Values: [tty, plain, html]
    --readName
      Show read name
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -r, --region
      Interval list
      Default: []
    --samFilter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -V, --variant, --variants, --vcf
      Variant file. if filename ends with '.list' it is interpreted as a list 
      of file (one file per line)
    --version
      print version and exit
    -width, --width
      default screen width
      Default: -1

```


## Keywords

 * sam
 * bam
 * visualization
 * terminal


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew tview
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/tview/TViewCmd.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/tview/TViewCmd.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/tview/TViewCmdTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/tview/TViewCmdTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **tview** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/tview.jar -R toy.fa toy.bam --clip  --groupby sample  --insert  --plain

ref:1-61
          POS: 1.......^^..11..^^^^....^^..21........31...^.....41........51
          REF: AGCATGTT**AGATAA****GATA**GCTGTGCTAGTAGGCAG*TCAGCGCCATNNNNNNN
               
ndefined_sample      TT**AGATAAAGAGGATA**-CTG               cagcgccat       
                      AAAAGATAAGG**GATAAA    NNNNNNtaggc                    
                  NNNNN**AGCTAA                                             
                                    ATA**GCT--------------CTCAGC            
ample CONSENSUS      TTAAAGNTAANGAGGATAAAGCTG      TAGGC  CTCAGCGCCAT
efined_sample 3          ********  ****   **                ****     
                         ********  ****   **                ****     
                         ********  ****   **                ****     
                     ************************      *****  ***********
                     ************************      *****  ***********
                     ************************      *****  ***********
               ******************************************************
               ******************************************************
               ******************************************************
             0 ******************************************************
```
            
