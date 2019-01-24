# VcfTail

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

print the last variants of a vcf


## Usage

```
Usage: vcftail [options] Files
  Options:
    -c, --bycontig
      number of variants
      Default: false
    -n, --count
      number of variants
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcftail
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfTail.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfTail.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfTailTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfTailTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcftail** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
java -jar dist/vcftail.jar -n2 |\
grep -v "##"| cut -f 1,2,4,5

#CHROM  POS REF ALT
chr1    935492  G   T
chr1    1334052 CTAGAG  C
```

