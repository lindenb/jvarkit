# SwingVcfView

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

VCFviewer using Java Swing UI


## Usage

```
Usage: swingvcfview [options] Files
  Options:
    --gff, --gff3
      GFF3 file Path used to find intervals by gene name
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -r, --regions, --interval
      default interval region on opening
      Default: <empty string>
    --limit
      Limit number of variants. Ignore if < 0
      Default: -1
    --pedigree
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    --version
      print version and exit

```


## Keywords

 * vcf
 * visualization
 * swing


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew swingvcfview
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210503

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingVcfView.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingVcfView.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingvcfview** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/swingvcfview.jar input.vcf.gz
```

## Screenshot

 * https://twitter.com/yokofakun/status/1392173413079322625

