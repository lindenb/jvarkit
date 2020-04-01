# ValidateCnv

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Experimental CNV Genotyping. Look variance of depths before/after putative known CNV.


## Usage

```
Usage: validatecnv [options] Files
  Options:
    -B, --bams, --bam
      Path to bam. File with suffix .list is interpretted as a file containing 
      a list of paths to bams.
      Default: []
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -x, --extend
      Search the boundaries in a region that is 'x'*(CNV-length). So if x if 
      0.5, a region chr1:100-200 will be searched chr1:50-100 + chr1:200-250
      Default: 0.5
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality.
      Default: 20
    --max, --max-size
      Max abs(SV) size to consider.A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 1000000
    --min-depth
      If minimum depth in region is lower than 'x', set the genotype as 
      NO_CALL 
      Default: 15
    --min-read-support-sv
      min number of read supporting SV.
      Default: 3
    --min, --min-size
      Min abs(SV) size to consider.A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
      Default: 50
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --stringency
      SAM Validation Stringency.
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    -t, --treshold
      HOM_DUP if 2.0-x<=depth<=2.0+x DUP if 1.5-x<=depth<=1.5+x . HET_DEL if 
      0.5-x<=depth<=0.5+x HOM_DEL if 0.0-x<=depth<=0.0+x . A decimal number 
      between 0.0 and 1.0. If the value ends with '%' it is interpretted as a 
      percentage eg. '1%' => '0.01'. A slash '/' is interpretted as a ratio. 
      e.g: '1/100' => '0.01'.
      Default: 0.05
    --version
      print version and exit

```


## Keywords

 * cnv
 * bam
 * sam
 * vcf
 * depth


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew validatecnv
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190914

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/ValidateCnv.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/ValidateCnv.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/ValidateCnvTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/ValidateCnvTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **validatecnv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a VCF (or a BED file

for VCF, only SVTYPE=DEL/INS/DUP are considered



## Example

```
find DIR -type f -name "*.bam" > bam.list
$ java -jar ${JVARKIT_HOME}/dist/validatecnv.jar -R reference.fa -B bam.list  20190320.MANTA.vcf 
```

