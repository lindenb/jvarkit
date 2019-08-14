# VcfGnomadSV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Peek annotations from gnomad structural variants


## Usage

```
Usage: vcfgnomadsv [options] Files
  Options:
    --any-overlap-filter
      If not empty, set this FILTER if any variant in gnomad is found 
      overlaping the variant BUT we didn't find a correct match
      Default: <empty string>
    --bnd-distance
      two bnd are identical if they're distant from  'x' bases
      Default: 10
    --discordant_svtype
      If not empty, set this FILTER if SVTYPE are discordants
      Default: <empty string>
    --fraction
      two segments are identical if they overlap at fraction 'x'
      Default: 0.7
  * -g, --gnomad
      Gnomad-SV VCF file. see 
      https://gnomad.broadinstitute.org/downloads#structural-variants 
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --in-gnomad-filter
      If not empty, set this FILTER is variant was found in gnomad
      Default: <empty string>
    -o, --output
      Output file. Optional . Default: stdout
    -p, --prefix
      INFO field prefix
      Default: GNOMAD_
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * gnomad
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgnomadsv
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190814

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSV.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSVTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gnomad/VcfGnomadSVTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgnomadsv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


