# NFHashDiff

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compare .nextflow.log


## Usage

```
Usage: java -jar dist/nfhashdiff.jar  [options] Files
Usage: nfhashdiff [options] Files
  Options:
    -g, --grep
      Search for all those strings
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -1
      hide unique to file1
      Default: false
    -2
      hide unique to file2
      Default: false
    -t
      look only at hashes value that value

```


## Keywords

 * json


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew nfhashdiff
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20251119

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/nextflow/hashdiff/NFHashDiff.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/nextflow/hashdiff/NFHashDiff.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **nfhashdiff** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Nextflow run **must** be launcher with `-dump-hashes json`

