# SamCustomSortJdk

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sort a BAM file using a java expression compiled at runtime.


## Usage

```
Usage: samcustomsortjdk [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
      java expression
    -f, --file
      java file. Either option -e or -f is required.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --saveCodeInDir
      Save the generated java code in the following directory
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * java
 * jdk
 * sort



## See also in Biostars

 * [https://www.biostars.org/p/305181](https://www.biostars.org/p/305181)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samcustomsortjdk
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamCustomSortJdk.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samjs/SamCustomSortJdk.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samjs/SamCustomSortJdkTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/samjs/SamCustomSortJdkTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samcustomsortjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Sort a BAM using a java expression compiled at runtime

## How it works

The user provides some java code of a `ComparatorSamRecord>`. This class extends `SamCustomSortJdk.AbstractSamComparator`.
At the time of writting, the  `SamCustomSortJdk.AbstractSamComparator` is:


```java
public static abstract class AbstractSamComparator
	implements Comparator<SAMRecord>
	{
	//input SAM header
	protected final SAMFileHeader header;
	//coordinate comparator 
	private final SAMRecordCoordinateComparator coordinateComparator;
	// query name comparator
	private final SAMRecordQueryNameComparator queryNameComparator;
	protected AbstractSamComparator(final SAMFileHeader header) {
		this.header = header;
		this.coordinateComparator = new SAMRecordCoordinateComparator();
		this.queryNameComparator = new SAMRecordQueryNameComparator();
		}
	// get a SAMRecordCoordinateComparator ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordCoordinateComparator.html ), never null 
	protected SAMRecordCoordinateComparator getCoordinateComparator() {
		return this.coordinateComparator;
	}
	// get a SAMRecordQueryNameComparator, never null 
	protected SAMRecordQueryNameComparator getQueryNameComparator() {
		return this.queryNameComparator;
	}
	// sort reads using 'SAMRecordCoordinateComparator.fileOrderCompare' https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordCoordinateComparator.html
	protected int fileOrderCompare(final SAMRecord R1, final SAMRecord R2) {
    	return getCoordinateComparator().fileOrderCompare(R1, R2);
    	}
	}
```
 
The user's code can be the code of the function `int compareTo(final SAMRecord R1,final SAMRecord R2)`, for example to sort on the mapping quality:

```java
return R1.getMappingQuality() - R2.getMappingQuality();
```

or the whole body of the comparator (option --body )


```java
private int mapq(final SAMRecord R) {
return R.getMappingQuality();
}

public int compareTo(final SAMRecord R1,final SAMRecord R2) {
return mapq(R1) - mapq(R2);
}
```

## Example


sort on amount of reference sequence covered, using the cigar string

```
 java -jar dist/samcustomsortjdk.jar  --body -e 'private int score(final SAMRecord R) { if(R.getReadUnmappedFlag() || R.getCigar()==null) return 0; return R.getCigar().getReferenceLength();} @Override public int compare(SAMRecord a,SAMRecord b) { return Integer.compare(score(a),score(b));}' in.bam
 ```

## History

 * 2019-01: migrating to openjkd: switched to in-memory compiling to external compiling

