# FastqRecordTreePack

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Create a TreeMap from one or more Fastq file. Ouput is a SVG file


## Usage

```
Usage: fastqrecordtreepack [options] Files
  Options:
    -c, --config
      XML config file
    -x, --dimension
      dimension of the output rectangle
      Default: 1000x1000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew fastqrecordtreepack
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/treepack/FastqRecordTreePack.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/treepack/FastqRecordTreePack.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqrecordtreepack** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Synopsis




```
$ java -jar dist/fastqrecordtreepack.jar -c config.xml (stdin|fq1.gz fq2.gz ...) gt; out.svg
```





### XML config


XML root is <treepack>. children is '<node>' .
A '<node>' has an attribute 'name'. The text content of the <node> will be evaluated as a javascript expression with the embedded javascript engine.
The javascript engine injects record a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html and
header a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html.



### Example



```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
	<node name="length">record.length()</node>
	<node name="firstBase">(record.length()&gt;0?record.getReadString().charAt(0):null)</node>
</treepack>

```

s


```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA21144/sequence_read/ERR047877.filt.fastq.gz" |\
   gunzip -c | java -jar dist/fastqrecordtreepack.jar -c) config.xml  > out.svg

```



![img](https://pbs.twimg.com/media/Bem-_tVCEAA9uT1.jpg:large)



### See also



 *  VcfTreePack
 *  http://www.cs.umd.edu/hcil/treemap-history/







### XML config

XML root is <treepack>. children is '<node>' .
A '<node>' has an attribute 'name'. The text content of the <node> will be evaluated as a javascript expression with the embedded javascript engine.
The javascript engine injects **record** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/fastq/FastqRecord.html) and
**header** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMFileHeader.html).


### Example

```
$ cat config.xml

<?xml version="1.0"?>
<treepack>
	<node name="length">record.length()</node>
	<node name="firstBase">(record.length()&gt;0?record.getReadString().charAt(0):null)</node>
</treepack>

```
s

```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA21144/sequence_read/ERR047877.filt.fastq.gz" |\
   gunzip -c | java -jar dist/fastqrecordtreepack.jar -c) config.xml  > out.svg

```


![https://pbs.twimg.com/media/Bem-_tVCEAA9uT1.jpg:large](https://pbs.twimg.com/media/Bem-_tVCEAA9uT1.jpg:large)



