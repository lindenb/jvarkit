# Biostar497922

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF into separate VCFs by SNP count


## Usage

```
Usage: biostar497922 [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -C, --chrom
      by chromosomes when using '-n'.
      Default: false
    --count, -n
      number of variants per vcf
      Default: -1
    -d, --distance
      max distance beween consecutive variants (ignore if <=0) . A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: -1
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --index
      start numering file index from 'x'
      Default: 1
    -D, --length
      max distance beween first and last variant (ignore if <=0) . A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: -1
    -m, --manifest
      output manifest to this file.
  * -o, --output
      Output directory
    --prefix
      file prefix
      Default: split
    --version
      print version and exit

```


## Keywords

 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/497922](https://www.biostars.org/p/497922)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar497922
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210319

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar497922.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar497922.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar497922** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
```


