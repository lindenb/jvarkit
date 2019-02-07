# LumpyVcfToCircos

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Lumpy to Circos


## Usage

```
Usage: lumpyvcf2circos [options] Files
  Options:
    -bnb, -bnd, --bnd
      hide SVTYPE=BND
      Default: false
    -del, --del, --deletions
      hide <DEL>
      Default: false
    -dup, --dup, --duplications
      hide <DUP>
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -inv, --inv, --invertions
      hide <INV>
      Default: false
    --su, --minsu
      Min supporting reads.
      Default: 20
    -o, --output
      output directory or zip file
    -p, --prefix
      file prefix
      Default: lumpy.
    --version
      print version and exit

```


## Keywords

 * lumpy
 * circos
 * sv
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew lumpyvcf2circos
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/lumpysv/LumpyVcfToCircos.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/lumpysv/LumpyVcfToCircos.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **lumpyvcf2circos** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
 java -jar dist/lumpyvcf2circos.jar  --minsu 50 -inv -bnb -dup  -o tmp  LumpyExpress.vcf.gz \
  && (cd tmp; /path/to/bin/circos  -outputdir . -conf lumpy.circos.conf  )
```

