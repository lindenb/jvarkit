# VcfPeekVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Get the INFO from a VCF and use it for another VCF


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfpeekvcf  [options] Files

Usage: vcfpeekvcf [options] Files
  Options:
    -a, -alt, --alt
      How alt allele must be found in the variants of the indexed file.All: 
      All ALT alleles must be found in the database ALTs. at_least_one: At 
      least one of the user ALT must be found in database ALTs. None: just use 
      CHROM/POS/REF 
      Default: none
      Possible Values: [none, all, at_least_one]
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -b, --buffer-size
      When we're looking for variants in a lare VCF file, load the variants in 
      an interval of 'N' bases instead of doing a random access for each 
      variant.. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 100000
    --default-float
      default value for Type=Float
    --default-int
      default value for Type=Integer
    --default-string
      default value for Type=String
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -missingIsError, --missingIsError
      Missing Info Header is an error
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -p, --prefix
      prefix all database tags with this prefix to avoid collisions
      Default: <empty string>
    -i, --replaceid
      Replace the ID field if it exists
      Default: false
  * -f, --tabix, --resource
      The VCF file indexed with TABIX or tribble. Source of the annotations
    -span, --span
      [20180713] when checking for the '--alt' option, ignore spanning 
      deletion: *
      Default: false
    -t, --tags
      tag1,tag2,tag... the INFO keys to peek from the indexed file
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation



## Creation Date

20150521

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpeekvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


