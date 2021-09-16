# RNASeqPolyA

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

find poly-A tail in RNASeq data


## Usage

```
Usage: rnaseqpolya [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -C, --contig
      limit to this contig/chromosome
    --disable-index
      Disable use of BAM index
      Default: false
    -d, --duplicate-ends
      keep only one transcript if transcripts share the same 3' coordinate
      Default: false
    --filter-reads
      remove duplicate, supplementary, malformed reads
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -g, --gff, --gff3
      GFF3 file containing the exons.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -indels, --indels
      Ignore reads containing indels.
      Default: false
    -o, -out, --out
      Output file. Optional . Default: stdout
    -x, --overlapping-exon
      Ignore exon if an exon from another transcript overlaps the end of the 
      last exon.
      Default: false
    -p, --primer
      Search for poly-A 'A{x}' dandling part of the read. -1 : ignore and 
      search for poly-A just after the exon boundary . If it's found, we count 
      the number of A starting with this pattern.
      Default: -1
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * bam
 * sam
 * rnaseq
 * polya


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew rnaseqpolya
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210913

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/rnaseqpolya/RNASeqPolyA.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/rnaseqpolya/RNASeqPolyA.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **rnaseqpolya** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

find evidence of poly-A tails in RNASeq data.

## Input

input is a set of path to indexed BAM/CRAM files or a file with the `.list` suffix containing the path to the BAM/CRAM files (one per line)

## Output

ouput is a VCF file. Each variant is a transcript.

## Example

```bash
$ find dir1 -type f -name "*.bam" > in.list
$ java -jar rnaseqpolya.jar -p 5 --reference ref.fa --gff3 in.gff  --out out.vcf.gz  in.list 
```

## See also

https://twitter.com/yokofakun/status/1438137720484900869

![twitter](https://pbs.twimg.com/media/E_VKkVJWEAMA6xy?format=png&name=900x900 "Screenshot")

