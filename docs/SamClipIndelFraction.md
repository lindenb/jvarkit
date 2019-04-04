# SamClipIndelFraction

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract clipping/indel fraction from BAM


## DEPRECATED

This tool can be replace with Bioalcidaejdk

## Usage

```
Usage: Samclipindelfraction [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -t
      type
      Default: allclip
      Possible Values: [leftclip, rightclip, allclip, insert, deletion]

```


## Keywords

 * sam
 * bam
 * clip


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew Samclipindelfraction
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFraction.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFraction.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFractionTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFractionTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **Samclipindelfraction** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ samtools view -h -F3844 my.bam  | java -jar dist/samclipindelfraction.jar 

##UNMAPPED_READS=0
##MAPPED_READS=3028359
##CLIPPED_READS=1182730
##CLIPPED_READS_5_PRIME=597757
##CLIPPED_READS_3_PRIME=617399
##UNCLIPPED_READS=1845629
##COUNT_BASES=338644685
#CLIP	COUNT	FRACTION_OF_MAPPED_READS
0	1845629	0.5
1	7	1.8963724562195327E-6
2	6756	0.0018302703306027376
3	695	1.8828269386751074E-4
4	794	2.1510281860547272E-4
5	819	2.2187557737768533E-4
6	471	1.275987752684857E-4
7	447	1.210969268471616E-4
(...)
```

plotting:
```bash
$ java -jar dist/samclipindelfraction.jar |\
   grep -v "##" | cut -f1,2 | tr -d '#' > output.txt
```

then, in R:
```R
T<-read.table('output.txt',header=TRUE)
plot(T[T$CLIP>0,])
```

