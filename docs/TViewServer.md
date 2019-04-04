# TViewServer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Web Server displaying SAM/BAM file. A web interface for jvarkit:tview


## Usage

```
Usage: tviewserver [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --max
      Max interval Length
      Default: 2000
    -nojs, --no-javascript
      Disable Javascript (which is not filesystem-safe).
      Default: false
    -P, --port, -port
      Server listening port
      Default: 8080
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --shutdown-after
      Stop the server after 'x' seconds.
      Default: -1
    --url
      A custom URL for a web browser. The following words will be replaced by 
      their values: ${CHROM}, ${START}, ${END}. For example for IGV that would 
      be: 'http://localhost:60151/goto?locus=${CHROM}%3A${START}-${END}' (see 
      http://software.broadinstitute.org/software/igv/book/export/html/189) 
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * table
 * visualization
 * server
 * web


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew tviewserver
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/tview/TViewServer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/tview/TViewServer.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/tview/TViewServerTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/tview/TViewServerTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **tviewserver** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

Input is a set of indexed BAM files  or a file containing the path to the BAMs.

## Screenshot

https://twitter.com/yokofakun/status/925395272225746944

![twitter](https://pbs.twimg.com/media/DNepyjOW4AA-_rz.jpg "Screenshot")


## Example 

```
$ java -jar dist/tviewserver.jar -R src/test/resources/toy.fa src/test/resources/toy.bam
2017-10-31 16:31:15.281:INFO::main: Logging initialized @1262ms
[INFO][TViewServer]Starting com.github.lindenb.jvarkit.tools.tview.TViewServer on http://localhost:8080
2017-10-31 16:31:15.523:INFO:oejs.Server:main: jetty-9.3.7.v20160115
2017-10-31 16:31:15.690:INFO:oejs.ServerConnector:main: Started ServerConnector@14a2189{HTTP/1.1,[http/1.1]}{0.0.0.0:8080}
2017-10-31 16:31:15.691:INFO:oejs.Server:main: Started @1675ms

```


