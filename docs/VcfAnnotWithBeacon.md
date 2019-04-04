# VcfAnnotWithBeacon

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with ga4gh beacon


## Usage

```
Usage: vcfannotwithbeacon [options] Files
  Options:
    --baseurl
      Beacon Base URL API
      Default: https://beacon-network.org/api
    -B, --bdb
      Optional BerkeleyDB directory to store result. Avoid to make the same 
      calls to beacon
    --build
      genome build
      Default: HG19
    --cert
      ignore SSL certification errors
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --noupdate
      Don't query the variant already having the tag / do not update the 
      existing annotation
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --stopOnError
      Stop on network error.
      Default: false
    --tag, -T
      INFO TAG
      Default: BEACON
    --tee
      show what's happening in the network
      Default: false
    --version
      print version and exit

```


## Keywords

 * ga4gh
 * beacon
 * vcf
 * annotation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfannotwithbeacon
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ga4gh/VcfAnnotWithBeacon.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/ga4gh/VcfAnnotWithBeacon.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/ga4gh/VcfAnnotWithBeaconTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/ga4gh/VcfAnnotWithBeaconTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfannotwithbeacon** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 
