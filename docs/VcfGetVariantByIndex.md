# VcfGetVariantByIndex

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Access a Plain or BGZF-compressed VCF file by index


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfbyindex  [options] Files

Usage: vcfbyindex [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -f, --force
      Force the creation of the index
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -I, --index
      Comma separated of 1-based index for query
      Default: <empty string>
    -i, --index-file
       (file) list of 1-based indexes for query
    --version
      print version and exit
    -o
      Output file. Optional . Default: stdout

```


## Keywords

 * vcf


## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbyindex/VcfGetVariantByIndex.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbyindex/VcfGetVariantByIndex.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbyindex/VcfGetVariantByIndexTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfbyindex/VcfGetVariantByIndexTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbyindex** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
# get random indexes
$ gunzip -c input.vcf.gz | grep -v "#" | awk '{print NR;}' | shuf | head -n 15 > index.list
# get those 15 variants
java -jar dist/jvarkit.jar vcfgetvariantbyindex -i index.list input.vcf.gz > output.vcf

```


