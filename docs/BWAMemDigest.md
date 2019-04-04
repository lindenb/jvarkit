# BWAMemDigest

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)




## Usage

```
Usage: bwamemdigest [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --version
      print version and exit
    -B
      BED of Regions to ignore.
    -x
       Extend BED Regions to ignore by 'x' bases.
      Default: 0

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bwamemdigest
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mem/BWAMemDigest.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mem/BWAMemDigest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bwamemdigest** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

package/bwa/bwa-0.7.4/bwa mem \
 -M  ~/tmp/DATASANGER/hg18/chr1.fa  \
 ~/tmp/DATASANGER/4368_1_1.fastq.gz  ~/tmp/DATASANGER/4368_1_2.fastq.gz 2> /dev/null | 
 java -jar dist/bwamemdigest.jar -B <(curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" | gunzip -c | cut -f2,3,4 ) -x 500  | tee /dev/tty | gzip --best > /tmp/jeter.mem.bed.gz


