# VcfRemoveUnusedAlt

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Remove unused ALT allele if there is no genotype with this alt, or there is no sample but AC=0


## Usage

```
Usage: vcfremoveunusedalt [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -neverspan, --neverspan
      Remove ALL spanning deletions '*'. VCF must have no genotype.
      Default: false
    -onespan, --onespan
      Don't print the variant if the only remaining allele is  '*'
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * genotype


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfremoveunusedalt
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveUnusedAlt.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveUnusedAlt.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveUnusedAltTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfRemoveUnusedAltTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfremoveunusedalt** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

when using gatk SelectVariants with sample names (-sn) some alleles specific of the samples than have been removed, remain in the vcf.

## SNPEFF / VEP

this tool removes unused annotations from SNPEFF(ANN=) and VEP.

## Example

```bash
$ cat in.vcf
(...)
chr1	7358	.	ACTT	*,A	1313.61	PASS	AC=0,10;AF=0,0.005828;AN=1716

$ java -jar dist/vcfremoveunusedalt.jar  in.vcf | grep -w 17358 -m1
chr1	7358	.	ACTT	A	1313.61	PASS	AC=10;AF=0.005828;AN=1716
```
