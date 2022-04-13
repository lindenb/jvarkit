# SwingVcfJexlFilter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter VCF using Java Swing UI and JEXL/Javascript expression


## Usage

```
Usage: java -jar dist/swingvcfjexl.jar  [options] Files
Usage: swingvcfjexl [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --version
      print version and exit

```


## Keywords

 * vcf
 * visualization
 * swing
 * jexl
 * javascript


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew swingvcfjexl
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20220413

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingVcfJexlFilter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingVcfJexlFilter.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingvcfjexl** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/swingvcfjexl.jar src/test/resources/rotavirus_rf.vcf.gz
```

## Screenshot

https://twitter.com/yokofakun/status/1514287386410336259

![https://twitter.com/yokofakun/status/1514287386410336259](https://pbs.twimg.com/media/FQPUMpbWUAgJbCs?format=png&name=900x900)

