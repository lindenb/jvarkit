# SortSamRefName

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Sort a BAM of contig and then on name


## Usage

```
Usage: sortsamrefname [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
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
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * sam
 * sort



## See also in Biostars

 * [https://www.biostars.org/p/154220](https://www.biostars.org/p/154220)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew sortsamrefname
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SortSamRefName.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SortSamRefName.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SortSamRefNameTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SortSamRefNameTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sortsamrefname** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$  java -jar dist/sortsamrefname.jar /commun/data/packages/samtools/1.2/samtools/examples/toy.sam  2> /dev/null 
@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@CO	SortSamRefName 1c7bc5e674136947586779a2aac53e576db4a67f /commun/data/packages/samtools/1.2/samtools/examples/toy.sam
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*
r003	0	ref	9	30	5H6M	*	0	0	AGCTAA	*
r003	16	ref	29	30	6H5M	*	0	0	TAGGC	*
r004	0	ref	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*
x1	0	ref2	1	30	20M	*	0	0	AGGTTTTATAAAACAAATAA	????????????????????
x2	0	ref2	2	30	21M	*	0	0	GGTTTTATAAAACAAATAATT	?????????????????????
x3	0	ref2	6	30	9M4I13M	*	0	0	TTATAAAACAAATAATTAAGTCTACA	??????????????????????????
x4	0	ref2	10	30	25M	*	0	0	CAAATAATTAAGTCTACAGAGCAAC	?????????????????????????
x5	0	ref2	12	30	24M	*	0	0	AATAATTAAGTCTACAGAGCAACT	????????????????????????
x6	0	ref2	14	30	23M	*	0	0	TAATTAAGTCTACAGAGCAACTA	???????????????????????
```

