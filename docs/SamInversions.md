# SamInversions

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan inversions in SAM using split and FF/RR reads on same strand.


## Usage

```
Usage: java -jar dist/saminversions.jar  [options] Files
Usage: saminversions [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -B, --bed
      Optional BED file. SV should overlap this bed.
    --control, --controls
      BAM/CRAM file that can be used as a control. I discordand reads are 
      found in any control, the candidate inversion is rejected. A file with 
      the '.list' suffix is interpreted as a file containing the path to the 
      bam files.
    --debug
      debug
      Default: false
    --distance
      search for breakend (BND) around 'x' bases.
      Default: 100
    --exclude-bed, --exclude
      skip reads pverlapping that bed file
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality.
      Default: 10
    --no-clip
      Ignore the clipping status of the SAMRecord, just look at supplementary 
      reads without looking at the clipping
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
      Default: []
    --version
      print version and exit
    --vcf-out, --out-vcf, --write-vcf
      save output in VCF format to that file

```


## Keywords

 * sam
 * bam
 * sv
 * inversion


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew saminversions
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200609

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/saminversions/SamInversions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/saminversions/SamInversions.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **saminversions** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




