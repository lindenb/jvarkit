# GtfRetroCopy

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan retrocopies by comparing the gtf/intron and the deletions in a VCF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gtfretrocopy  [options] Files

Usage: gtfretrocopy [options] Files
  Options:
    --all
      all introns must be found
      Default: false
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -d, --distance
      max distance between an intron and the deletion found in the VCF
      Default: 10
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -gtf, --gtf
      A GTF (General Transfer Format) file. See 
      https://www.ensembl.org/info/website/upload/gff.html . Please note that 
      CDS are only detected if a start and stop codons are defined.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --id, -id
      Which key should I use for the column ID. The idea is to use the gene 
      name to get the uniq entities per vcf.
      Default: transcript_id
      Possible Values: [transcript_id, gene_id, gene_name]
    -k, --known
      Gene-ID of known retrogenes. One per line. A source could be : 
      http://retrogenedb.amu.edu.pl/static/download/ 
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * gtf
 * retrocopy
 * deletion



## Creation Date

20190813

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/retrocopy/GtfRetroCopy.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/retrocopy/GtfRetroCopy.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/retrocopy/GtfRetroCopyTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/retrocopy/GtfRetroCopyTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gtfretrocopy** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/gtfretrocopy.jar --gtf transcript.gtf.gz input.vcf.gz > retrocopies.vcf
```


