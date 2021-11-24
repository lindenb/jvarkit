# VcfROH

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

VCF ROH


## Usage

```
Usage: vcfroh [options] Files
  Options:
    --bed
      Output BED instead of interval list
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --min-length
      Minimum block length. A distance specified as a positive integer.Commas 
      are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 1000
    --min-score
      Minimum block score.
      Default: 0.0
    --min-variants
      Minimum number of variants.
      Default: 10
    --nc2hr
      treat NO_CALL to HOM_REF
      Default: false
    --output, -o
      Output file. Optional . Default: stdout
    --score
      HOM_REF/HOM_VAR score. HET score will be 'x' - 1.0.
      Default: 0.0025
    --version
      print version and exit

```


## Keywords

 * vcf
 * roh


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfroh
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20211123

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/roh/VcfROH.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/roh/VcfROH.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfroh** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


 
## Example


```bash
java -jar dist/vcfpar.jar   in.vcf > out.vcf
```

