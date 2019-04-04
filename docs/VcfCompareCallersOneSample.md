# VcfCompareCallersOneSample

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

For my colleague Julien: VCF with one sample called using different callers. Only keep variant if it was found in min<x=other-files<=max


## Usage

```
Usage: vcfcomparecallersonesample [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -M
       max number of challengers found, inclusive.
      Default: 2147483646
    -a
      ignore ALT allele
      Default: false
    -f
      VCF to be challenged.  Must be sorted on dict. Must contain a dict.
      Default: []
    -m
      min number of challengers found, inclusive.
      Default: 0

```


## Keywords

 * vcf
 * compare


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcomparecallersonesample
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallersOneSample.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallersOneSample.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomparecallersonesample** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Misc.

I used 
```
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f samtools.vcf gatk.vcf
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f gatk.vcf samtools.vcf
```

shouldn't I get the same number of variants in both files ?
**Answer** is "not always" in the following case:

in gatk.vcf:

```
11	244197	rs1128322	T	C
````

in samtools.vcf:
```
11	244197	rs1128322	T	C,G
```

in
```
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f samtools.vcf gatk.vcf
```
we keep the variant because we found 'C' in samtools and 'gatk'

```
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f gatk.vcf samtools.vcf
```

the variant is discarded because 'G' is found in samtools but not in 'gatk.vcf'

