# VcfGnomad

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Peek annotations from gnomad


## Usage

```
Usage: vcfgnomad [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --bufferSize
      When we're looking for variant in a lare VCF file, load the variants in 
      an interval of 'N' bases instead of doing a random access for each 
      variant. A distance specified as a positive integer.Commas are removed. 
      The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 10000
    -F, --fields
      AF fields to peek-up from gnomad. Space/comma/semicolon separated
      Default: AF_popmax,AF_nfe
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
  * -g, --gnomad
      Path to Indexed Gnomad VCF file.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --max-af
      Max allele frequency
      Default: 1.0
    --min-af
      Min allele frequency
      Default: 0.0
    --noUpdateId
      do Not Update ID if it is missing in user's variant
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --prefix
      If not empty, include the Gnomad FILTERs using this prefix.
      Default: GNOMAD
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * gnomad


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgnomad
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20170407

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomad.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomad.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomad** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
The tool was redesigned on July 2nd, 2020.

## Example:

```
java -jar dist/vcfgnomad.jar -g src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz src/test/resources/test_vcf01.vcf
```

