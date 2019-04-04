# VcfSkat

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Calculate SKAT score for a VCF.


## Usage

```
Usage: vcfskat [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -ped, --pedigree
      A pedigree is a text file delimited with tabs. No header. Columns are 
      (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' 
      (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown  If not defined, I will try to extract the pedigree 
      from  the VCFheader.
    -p, --print
      Don't print the VCF on output,Just print the score on ouput, don't print 
      the VCF itself
      Default: false
    --skat-accept-filtered
      accept variants FILTER-ed
      Default: false
    --skat-adjusted
      SKAT adjusted
      Default: false
    --skat-num-retry
      compute n-times the p-value
      Default: 1
    --skat-optimized
      SKAT optimized (SKATO)/ davies method.
      Default: false
    --skat-random-seed
      Rstats value for `set.seed`. -1 == use random
      Default: -1
    -tag, --tag
      tag name in VCF header "tag=p-value"
      Default: SKAT
    --version
      print version and exit

```


## Keywords

 * vcf
 * pedigree
 * skat
 * burden


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfskat
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/skat/VcfSkat.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/skat/VcfSkat.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfskat** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/vcfskat.jar vcf_with_samples.vcf | grep SKAT
##SKAT=0.2215079
```

```
$ java -jar dist/vcfskat.jar -p vcf_with_samples.vcf 
0.2215079
```


