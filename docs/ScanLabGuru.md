# ScanLabGuru

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

scan the files stored in labguru


## Usage

```
Usage: java -jar dist/scanlabguru.jar  [options] Files
Usage: scanlabguru [options] Files
  Options:
  * --api-key, -A
      API key (can be a string or an existing file
      Default: <empty string>
    --base
      Labguru base URI.
      Default: https://cle.inserm.fr
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * labguru
 * lims
 * vcf
 * sam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew scanlabguru
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20240325

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/labguru/ScanLabGuru.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/labguru/ScanLabGuru.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **scanlabguru** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Motivation: scan the files stored in labguru


usage:

```
jvarkit -jar jvarkit.jar scanlabguru -A api.txt <projectid1> <projectid2> <projectid3>
```




